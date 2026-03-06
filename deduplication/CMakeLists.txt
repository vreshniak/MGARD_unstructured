cmake_minimum_required(VERSION 3.5)

#==================================================================================================
project(node_deduplication CXX)

find_package(PkgConfig REQUIRED)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
# pkg_search_module(MPI_CXX REQUIRED mpi-cxx)

# check for ADIOS
find_package(ADIOS2 PATHS /ccs/home/reshniakv/libs_frontier/adios2/lib64 REQUIRED COMPONENTS CXX)
# find_package(ADIOS2 REQUIRED PATHS /home/qian/Software/ADIOS2/install-adios/)
# find_package(ADIOS2 PATHS /ccs/home/reshniakv/libs_andes/adios/lib64 REQUIRED COMPONENTS CXX)



# Collect list of all source files
#==================================================================================================
# E.g: find src -type f \( -name "*.cc" -o -name "*.h" -o -name "*.hpp" \) | sort
set(SRC_DEDUP_FILES ./dedup.cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
# Sets the output directory for all executables to "bin"
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(dedup ${SRC_DEDUP_FILES})
# target_include_directories(dedup PRIVATE include)
target_include_directories(dedup PRIVATE .)
target_link_libraries(dedup PRIVATE adios2::cxx_mpi MPI::MPI_CXX OpenMP::OpenMP_CXX)
