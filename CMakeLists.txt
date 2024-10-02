cmake_minimum_required(VERSION 3.5)

#==================================================================================================
project(unstructured_serialization CXX)

find_package(MPI REQUIRED)

# check for ADIOS
find_package(ADIOS2 PATHS /ccs/home/reshniakv/libs_frontier/adios/lib64 REQUIRED COMPONENTS CXX)
# find_package(ADIOS2 REQUIRED PATHS /home/qian/Software/ADIOS2/install-adios/)
# find_package(ADIOS2 PATHS /ccs/home/reshniakv/libs_andes/adios/lib64 REQUIRED COMPONENTS CXX)

# check for MGARD
find_package(mgard PATHS /ccs/home/reshniakv/MGARD/install-hip-frontier REQUIRED COMPONENTS CXX)


# Collect list of all source files
#==================================================================================================
# E.g: find src -type f \( -name "*.cc" -o -name "*.h" -o -name "*.hpp" \) | sort
set(SRC_DATA_FILES
	src/extract_data.cpp
)
set(SRC_MERGE_FILES
	src/merge.cpp
)
set(SRC_REORDER_FILES
	src/reorder.cpp
)
set(SRC_COMPRESS_MGARD_FILES
	src/compress_mgard.cpp
)
set(SRC_COMPRESS_ADIOS_FILES
	src/compress_adios.cpp
)

add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(extract_data ${SRC_DATA_FILES})
target_include_directories(extract_data PRIVATE include)
target_link_libraries(extract_data PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(merge ${SRC_MERGE_FILES})
target_include_directories(merge PRIVATE include)
target_link_libraries(merge PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(reorder ${SRC_REORDER_FILES})
target_include_directories(reorder PRIVATE include)
target_link_libraries(reorder PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(compress_mgard ${SRC_COMPRESS_MGARD_FILES})
target_include_directories(compress_mgard PRIVATE include)
target_link_libraries(compress_mgard PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(compress_adios ${SRC_COMPRESS_ADIOS_FILES})
target_include_directories(compress_adios PRIVATE include)
target_link_libraries(compress_adios PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)