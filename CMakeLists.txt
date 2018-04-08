cmake_minimum_required (VERSION 3.1)

# Project name
project (quest)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
            "Choose the type of build, options are: None Debug Release."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

option(USE_CKB "Use checkboard decomposition" OFF)
#option(USE_CUDA "Use GPU support" OFF)

enable_language(Fortran)

# set cmake module path for *.cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

add_library(util OBJECT SRC/dqmc_mpi.F90 SRC/dqmc_util.F90 SRC/dqmc_cfg.F90 SRC/dqmc_omp.F90)
add_library(latt OBJECT SRC/dqmc_geom_param.F90 SRC/dqmc_latt.F90 SRC/dqmc_bonds.F90 SRC/dqmc_reclatt.F90 )
add_library(hamilt OBJECT SRC/dqmc_hamilt.F90 SRC/dqmc_symm.F90)
if (USE_CKB)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDQMC_CKB")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDQMC_CKB")
    add_library(mat OBJECT SRC/dqmc_wspace.F90 SRC/dqmc_struct.F90 SRC/dqmc_geom_wrap.F90 SRC/dqmc_checkerboard.F90)
else()
    add_library(mat OBJECT SRC/dqmc_wspace.F90 SRC/dqmc_struct.F90 SRC/dqmc_geom_wrap.F90 SRC/dqmc_matb.F90)
endif()
add_library(phy OBJECT SRC/dqmc_seqb.F90 SRC/dqmc_gfun.F90 SRC/dqmc_phy0.F90 SRC/computeg.c++
        SRC/blaslapack.c++
        SRC/kernel_cpu.c++
        SRC/profile.c++)
add_library(phy2 OBJECT SRC/dqmc_phy2.F90 SRC/dqmc_hubbard.F90 SRC/dqmc_gtau.F90)
add_library(tdm OBJECT SRC/dqmc_tdm1.F90 SRC/dqmc_tdm2.F90 SRC/dqmc_kbonds.F90 SRC/dqmc_2dperl.F90)

add_dependencies(latt util)
add_dependencies(hamilt latt)
add_dependencies(mat hamilt)
add_dependencies(phy mat)
add_dependencies(phy2 phy)
add_dependencies(tdm phy2)

add_library(dqmc STATIC  $<TARGET_OBJECTS:util>
                         $<TARGET_OBJECTS:latt>
                         $<TARGET_OBJECTS:hamilt>
                         $<TARGET_OBJECTS:mat>
                         $<TARGET_OBJECTS:phy>
                         $<TARGET_OBJECTS:phy2>
                         $<TARGET_OBJECTS:tdm>)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(ggeom EXAMPLE/geom/ggeom.F90)
target_link_libraries(ggeom dqmc)
target_include_directories(ggeom PUBLIC $<TARGET_PROPERTY:dqmc,Fortran_MODULE_DIRECTORY>)
set_target_properties(ggeom PROPERTIES LINKER_LANGUAGE Fortran)


if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -funroll-loops -Wall -O3 -fPIC")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -warn all -O3 -unroll -fPIC")
endif()

include_directories(${CMAKE_SOURCE_DIR}/SRC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -Wall -O3 -funroll-loops -fPIC")

find_package(MKL)
if("${MKL_LIBRARIES}" STREQUAL "")
    find_package(BLAS REQUIRED)
    if (BLAS_FOUND)
        message("MKL is not found using blas instead")
        target_link_libraries(dqmc blas)
    else()
        message(STATUS "No blas library could be found.")
    endif()
else()
    message(STATUS "Found required mkl libraries")
    target_link_libraries(dqmc ${MKL_LIBRARIES} pthread dl)
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
else()
    message(STATUS "OpenMP could not be found.")
endif()

install(FILES ${CMAKE_SOURCE_DIR}/src/py/MANIFEST.in DESTINATION ${CMAKE_BINARY_DIR})
install(FILES ${CMAKE_BINARY_DIR}/dqmc.so DESTINATION ${CMAKE_BINARY_DIR}/${PYDQMC_MOD})
install(FILES ${CMAKE_SOURCE_DIR}/src/py/dqmc/__init__.py DESTINATION ${CMAKE_BINARY_DIR}/${PYDQMC_MOD})
install(FILES ${CMAKE_SOURCE_DIR}/src/py/dqmc/dqmc_handler.py DESTINATION ${CMAKE_BINARY_DIR}/${PYDQMC_MOD})
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} setup.py install)")
