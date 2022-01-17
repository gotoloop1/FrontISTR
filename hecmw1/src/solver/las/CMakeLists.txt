###############################################################################
# Copyright (c) 2019 FrontISTR Commons
# This software is released under the MIT License, see License.txt
###############################################################################
cmake_minimum_required(VERSION 2.8.11)

list(APPEND hecmw_SOURCES
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_11.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_22.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_las.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_jadm.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_33.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_44.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_66.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling_nn.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_solver_scaling.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_tuning_fx.f90
  ${CMAKE_CURRENT_LIST_DIR}/crhcm.c
  ${CMAKE_CURRENT_LIST_DIR}/ebcm.c
  ${CMAKE_CURRENT_LIST_DIR}/first_element_sort.c
  ${CMAKE_CURRENT_LIST_DIR}/hcm.c
  ${CMAKE_CURRENT_LIST_DIR}/heap.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_reordering_helper.c
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_reordering_sort.f90
  ${CMAKE_CURRENT_LIST_DIR}/hecmw_reordering.f90
  ${CMAKE_CURRENT_LIST_DIR}/memory.c
  ${CMAKE_CURRENT_LIST_DIR}/nbcm.c
  ${CMAKE_CURRENT_LIST_DIR}/random.c
  ${CMAKE_CURRENT_LIST_DIR}/range_sum.c
  ${CMAKE_CURRENT_LIST_DIR}/rcm.c
  ${CMAKE_CURRENT_LIST_DIR}/single_sort.c
)

list(APPEND hecmw_INCLUDE_DIRS
  ${CMAKE_CURRENT_LIST_DIR}
)
