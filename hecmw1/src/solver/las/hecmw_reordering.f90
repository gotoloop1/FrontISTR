
module hecmw_reordering
    use hecmw_util
    use hecmw_matrix_misc
    use hecmw_reordering_sort

    implicit none

    private

    public :: hecmw_reordering_fw
    public :: hecmw_reordering_bw

    integer(kind=kint), allocatable :: order(:)
    integer(kind=kint), allocatable :: order_inverse(:)

    integer(kind=kint) :: n
    integer(kind=kint), allocatable :: edge_separator(:)
    integer(kind=kint), allocatable :: edge(:)
    integer(kind=kint), allocatable :: edge_weight(:)

contains

    subroutine hecmw_reordering_fw(hecMESH, hecMAT)
        type(hecmwST_local_mesh) :: hecMESH
        type(hecmwST_matrix) :: hecMAT

        integer(kind=kint) :: param1, param2

        write (*,*) 'hecmw_reordering_fw called'
        write (*,*) 'type:  ', hecmw_mat_get_reordering(hecMAT)

        if (hecmw_mat_get_reordering(hecMAT) /= 0) then
            call hecmw_construct_graph(hecMAT)
            
            allocate(order(n), order_inverse(n))
            
            select case (hecmw_mat_get_reordering(hecMAT))
            case (1)
                call reverse_cuthill_mckee_wrap(n, edge_separator, edge, edge_weight, order)
            case (2)
                call edge_based_cost_minimization_wrap(n, edge_separator, edge, edge_weight, order)
            case (3)
                call node_based_cost_minimization_wrap(n, edge_separator, edge, edge_weight, order)
            case (4)
                call hybrid_cost_minimization_wrap(n, edge_separator, edge, edge_weight, order)
            case (5)
                call coarsen_refine_hybrid_cost_minimization_wrap(n, edge_separator, edge, edge_weight, order)
            end select
            
            call hecmw_modify_reordering(hecMAT)

            call hecmw_check_order
            
            call hecmw_calculate_inverse
            
            call hecmw_destruct_graph
            
            call hecmw_apply_reordering(hecMESH, hecMAT)
        end if
    end subroutine hecmw_reordering_fw

    subroutine hecmw_reordering_bw(hecMESH, hecMAT)
        type(hecmwST_local_mesh) :: hecMESH
        type(hecmwST_matrix) :: hecMAT

        write (*,*) 'hecmw_reordering_bw called'

        if (hecmw_mat_get_reordering(hecMAT) /= 0) then
            call hecmw_restore_reordering(hecMESH, hecMAT)

            deallocate(order, order_inverse)
        end if
    end subroutine hecmw_reordering_bw

    subroutine hecmw_construct_graph(hecMAT)
        type(hecmwST_matrix) :: hecMAT

        integer(kind=kint), allocatable :: count(:)
        integer(kind=kint), allocatable :: buffer_separator(:)
        integer(kind=kint), allocatable :: buffer(:)
        integer(kind=kint) :: i, j, b, l, p

        n = hecMAT%NP

        allocate(edge_separator(n + 1))
        allocate(edge(2 * (hecMAT%indexL(n) + hecMAT%indexU(n))))
        allocate(edge_weight(2 * (hecMAT%indexL(n) + hecMAT%indexU(n))))
        allocate(count(n))
        allocate(buffer_separator(n + 1))
        allocate(buffer(2 * (hecMAT%indexL(n) + hecMAT%indexU(n))))

        count(:) = 0
        do i = 1, n
            b = hecMAT%indexL(i - 1) + 1
            l = hecMAT%indexL(i)
            count(i) = count(i) + (l - b + 1)
            do j = b, l
                p = hecMAT%itemL(j)
                count(p) = count(p) + 1
            end do
        end do
        do i = 1, n
            b = hecMAT%indexU(i - 1) + 1
            l = hecMAT%indexU(i)
            count(i) = count(i) + (l - b + 1)
            do j = b, l
                p = hecMAT%itemU(j)
                count(p) = count(p) + 1
            end do
        end do

        buffer_separator(1) = 0
        do i = 1, n
            buffer_separator(i + 1) = buffer_separator(i) + count(i)
            count(i) = 0
        end do

        do i = 1, n
            b = hecMAT%indexL(i - 1) + 1
            l = hecMAT%indexL(i)
            do j = b, l
                p = hecMAT%itemL(j)
                count(i) = count(i) + 1
                buffer(buffer_separator(i) + count(i)) = p - 1
                count(p) = count(p) + 1
                buffer(buffer_separator(p) + count(p)) = i - 1
            end do
        end do
        do i = 1, n
            b = hecMAT%indexU(i - 1) + 1
            l = hecMAT%indexU(i)
            do j = b, l
                p = hecMAT%itemU(j)
                count(i) = count(i) + 1
                buffer(buffer_separator(i) + count(i)) = p - 1
                count(p) = count(p) + 1
                buffer(buffer_separator(p) + count(p)) = i - 1
            end do
        end do

        edge_separator(1) = 0
        do i = 1, n
            b = buffer_separator(i) + 1
            l = buffer_separator(i + 1)
            edge_separator(i + 1) = edge_separator(i)
            call single_sort_wrap(l - b + 1, buffer(b:l), n - 1)
            p = -1
            do j = b, l
                if (p == buffer(j)) then
                    edge_weight(edge_separator(i + 1)) = edge_weight(edge_separator(i + 1)) + 1
                else
                    edge_separator(i + 1) = edge_separator(i + 1) + 1
                    edge(edge_separator(i + 1)) = buffer(j)
                    edge_weight(edge_separator(i + 1)) = 1
                    p = buffer(j)
                end if
            end do
        end do

        deallocate(count, buffer_separator, buffer)
    end subroutine hecmw_construct_graph

    subroutine hecmw_destruct_graph()
        deallocate(edge_separator, edge, edge_weight)
    end subroutine hecmw_destruct_graph

    subroutine hecmw_modify_reordering(hecMAT)
        type(hecmwST_matrix) :: hecMAT

        integer(kind=kint) :: tmp(hecMAT%NP - hecMAT%N)
        integer(kind=kint) :: i, count

        count = 0
        do i = 1, hecMAT%NP
            if (order(i) < hecMAT%N) then
                count = count + 1
                order(count) = order(i)
            else
                tmp(i - count) = order(i)
            end if
        end do
        order((hecMAT%N+1):hecMAT%NP) = tmp
    end subroutine hecmw_modify_reordering
    
    subroutine hecmw_calculate_inverse()
        integer(kind=kint) :: i

        do i = 1, n
            order_inverse(order(i) + 1) = i - 1
        end do
    end subroutine hecmw_calculate_inverse

    subroutine hecmw_apply_reordering(hecMESH, hecMAT)
        type(hecmwST_local_mesh) :: hecMESH
        type(hecmwST_matrix) :: hecMAT

        call hecmw_reordering_sort_system(hecMESH, hecMAT, order, order_inverse)
    end subroutine hecmw_apply_reordering

    subroutine hecmw_restore_reordering(hecMESH, hecMAT)
        type(hecmwST_local_mesh) :: hecMESH
        type(hecmwST_matrix) :: hecMAT

        call hecmw_reordering_sort_system(hecMESH, hecMAT, order_inverse, order)
    end subroutine hecmw_restore_reordering

    subroutine hecmw_check_order()
        integer(kind=kint) :: count(n)
        integer(kind=kint) :: i

        count = 0
        do i = 1, n
            count(order(i) + 1) = count(order(i) + 1) + 1
        end do

        do i = 1, n
            if (count(i) /= 1) then
                write(*,*) 'invalid order'
            end if
        end do
    end subroutine hecmw_check_order

end module hecmw_reordering
