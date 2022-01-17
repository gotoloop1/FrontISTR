
module hecmw_reordering_sort
    use hecmw_util
    use hecmw_matrix_misc

    implicit none

    private

    public hecmw_reordering_sort_system

contains

    subroutine hecmw_reordering_sort_system(hecMESH, hecMAT, order, order_inverse)
        type(hecmwST_local_mesh) :: hecMESH
        type(hecmwST_matrix) :: hecMAT
        integer(kind=kint) :: order(:)
        integer(kind=kint) :: order_inverse(:)

        call hecmw_reordering_sort_vector(hecMAT%NP, hecMAT%NDOF, order, hecMAT%X)
        call hecmw_reordering_sort_vector(hecMAT%NP, hecMAT%NDOF, order, hecMAT%B)

        call hecmw_reordering_sort_vector(hecMAT%NP, hecMAT%NDOF * hecMAT%NDOF, order, hecMAT%D)
        call hecmw_reordering_sort_off_diagonal(hecMAT, order, order_inverse)

        if (hecMESH%n_neighbor_pe /= 0) then
            call hecmw_reordering_sort_map_array(hecMESH%import_index(hecMESH%n_neighbor_pe), order_inverse, hecMESH%import_item)
            call hecmw_reordering_sort_map_array(hecMESH%export_index(hecMESH%n_neighbor_pe), order_inverse, hecMESH%export_item)
        end if
    end subroutine hecmw_reordering_sort_system

    subroutine hecmw_reordering_sort_vector(n, l, order, x)
        integer(kind=kint) :: n
        integer(kind=kint) :: l
        integer(kind=kint) :: order(:)
        real(kind=kreal) :: x(:)

        real(kind=kreal) :: buffer(n * l)
        integer(kind=kint) :: i, new_b, new_l, old_b, old_l

        do i = 1, n
            new_b = (i - 1) * l + 1
            new_l = i * l
            old_b = order(i) * l + 1
            old_l = (order(i) + 1) * l
            buffer(new_b:new_l) = x(old_b:old_l)
        end do
        x = buffer
    end subroutine hecmw_reordering_sort_vector

    subroutine hecmw_reordering_sort_off_diagonal(hecMAT, order, order_inverse)
        type(hecmwST_matrix) :: hecMAT
        integer(kind=kint) :: order(:)
        integer(kind=kint) :: order_inverse(:)

        integer(kind=kint), allocatable :: item_buffer_index(:)
        integer(kind=kint), allocatable :: item_buffer(:)
        real(kind=kreal), allocatable :: a_buffer(:)
        integer(kind=kint), allocatable :: tmp(:)
        integer(kind=kint) :: n, ndof2, i, j, k, old, b, l, p, count, l_count, u_count

        n = hecMAT%N
        ndof2 = hecMAT%NDOF * hecMAT%NDOF

        allocate(item_buffer_index(0:n))
        allocate(item_buffer(hecMAT%indexL(n) + hecMAT%indexU(n)))
        allocate(a_buffer((hecMAT%indexL(n) + hecMAT%indexU(n)) * ndof2))
        allocate(tmp(n * 2))

        item_buffer_index(0) = 0
        count = 0
        l_count = 0
        u_count = 0
        do i = 1, n
            old = order(i) + 1

            b = hecMAT%indexL(old - 1) + 1
            l = hecMAT%indexL(old)
            do j = b, l
                do k = 1, ndof2
                    a_buffer(count * ndof2 + k) = hecMAT%AL((j - 1) * ndof2 + k)
                end do
                p = order_inverse(hecMAT%itemL(j))
                item_buffer(count + 1) = p
                count = count + 1
                if (p + 1 < i) then
                    l_count = l_count + 1
                else
                    u_count = u_count + 1
                end if
            end do

            b = hecMAT%indexU(old - 1) + 1
            l = hecMAT%indexU(old)
            do j = b, l
                do k = 1, ndof2
                    a_buffer(count * ndof2 + k) = hecMAT%AU((j - 1) * ndof2 + k)
                end do
                p = order_inverse(hecMAT%itemU(j))
                item_buffer(count + 1) = p
                count = count + 1
                if (p + 1 < i) then
                    l_count = l_count + 1
                else
                    u_count = u_count + 1
                end if
            end do
            item_buffer_index(i) = count
        end do

        deallocate(hecMAT%itemL, hecMAT%itemU, hecMAT%AL, hecMAT%AU)
        allocate(hecMAT%itemL(l_count))
        allocate(hecMAT%itemU(u_count))
        allocate(hecMAT%AL(l_count * ndof2))
        allocate(hecMAT%AU(u_count * ndof2))

        hecMAT%indexL(0) = 0
        hecMAT%indexU(0) = 0
        do i = 1, n
            hecMAT%indexL(i) = hecMAT%indexL(i - 1)
            hecMAT%indexU(i) = hecMAT%indexU(i - 1)
            b = item_buffer_index(i - 1) + 1
            l = item_buffer_index(i)

            count = 0
            do j = b, l
                tmp(count * 2 + 1) = item_buffer(j)
                tmp(count * 2 + 2) = count
                count = count + 1
            end do
            call first_element_sort_wrap(l - b + 1, tmp, n - 1)

            count = 0
            do j = b, l
                if (tmp(count * 2 + 1) + 1 < i) then
                    old = b + tmp(count * 2 + 2) - 1
                    do k = 1, ndof2
                        hecMAT%AL(hecMAT%indexL(i) * ndof2 + k) = a_buffer(old * ndof2 + k)
                    end do
                    hecMAT%indexL(i) = hecMAT%indexL(i) + 1
                    hecMAT%itemL(hecMAT%indexL(i)) = tmp(count * 2 + 1) + 1
                else
                    old = b + tmp(count * 2 + 2) - 1
                    do k = 1, ndof2
                        hecMAT%AU(hecMAT%indexU(i) * ndof2 + k) = a_buffer(old * ndof2 + k)
                    end do
                    hecMAT%indexU(i) = hecMAT%indexU(i) + 1
                    hecMAT%itemU(hecMAT%indexU(i)) = tmp(count * 2 + 1) + 1
                end if
                count = count + 1
            end do
        end do

        deallocate(item_buffer, item_buffer_index, a_buffer, tmp)
    end subroutine hecmw_reordering_sort_off_diagonal

    subroutine hecmw_reordering_sort_map_array(n, order_inverse, x)
        integer(kind=kint) :: n
        integer(kind=kint) :: order_inverse(:)
        integer(kind=kint) :: x(:)

        integer(kind=kint) :: i

        do i = 1, n
            x(i) = order_inverse(x(i)) + 1
        end do
    end subroutine hecmw_reordering_sort_map_array

end module hecmw_reordering_sort
