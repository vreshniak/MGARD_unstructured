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

# find_package(OpenMP)

# find_package(zstd PATHS /ccs/home/reshniakv/MGARD/install-hip-frontier/lib/cmake/zstd/ REQUIRED)


find_package(PkgConfig)
# pkg_check_modules(ZSTD REQUIRED IMPORTED_TARGET GLOBAL zstd)
pkg_search_module(ZSTD IMPORTED_TARGET libzstd)

find_package(zfp PATHS /ccs/home/reshniakv/libs_frontier/zfp/lib64)

find_package(SZ3 REQUIRED PATHS /ccs/home/reshniakv/libs_frontier/sz3)
# find_package(SZ3 REQUIRED PATHS /lustre/orion/csc143/scratch/gongq/frontier/SoftwareDev/SZ3/install/usr/local/lib64/cmake/SZ3/)
# set(SZ3_INCLUDES  "/ccs/home/reshniakv/libs_frontier/sz3/include")
# set(ZSTD_INCLUDES "/ccs/home/reshniakv/MGARD/install-hip-frontier/include/")





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
set(SRC_ADD_REORDER_FILES
	src/add_ordering.cpp
)
set(SRC_MAKE_REORDER_FILES
	src/make_meshfile_with_ordering.cpp
)
set(SRC_COMPRESS_MGARD_FILES
	src/compress_mgard.cpp
)
set(SRC_COMPRESS_MGARD_SORT_FILES
	src/compress_mgard_sort.cpp
)
set(SRC_COMPRESS_ADIOS_FILES
	src/compress_adios.cpp
)
set(SRC_COMPRESS_ADIOS_PLUGIN_FILES
	src/compress_with_adios_plugin.cpp
)
set(SRC_DECOMPRESS_ADIOS_PLUGIN_FILES
	src/decompress_with_adios_plugin.cpp
)
set(SRC_EVALUATE_ERROR_FILES
	src/evaluate_compression_error.cpp
)
set(SRC_COMPRESS_ZFP_FILES
	src/compress_zfp.cpp
)
set(SRC_COMPRESS_SZ_FILES
	src/compress_sz.cpp
)


add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(extract_data ${SRC_DATA_FILES})
target_include_directories(extract_data PRIVATE include src)
target_link_libraries(extract_data PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(merge ${SRC_MERGE_FILES})
target_include_directories(merge PRIVATE include)
target_link_libraries(merge PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(reorder ${SRC_REORDER_FILES})
target_include_directories(reorder PRIVATE include)
target_link_libraries(reorder PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(add_ordering ${SRC_ADD_REORDER_FILES})
target_include_directories(add_ordering PRIVATE include)
target_link_libraries(add_ordering PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(make_meshfile_with_ordering ${SRC_MAKE_REORDER_FILES})
target_include_directories(make_meshfile_with_ordering PRIVATE include)
target_link_libraries(make_meshfile_with_ordering PRIVATE adios2::cxx11_mpi MPI::MPI_CXX)

add_executable(compress_mgard ${SRC_COMPRESS_MGARD_FILES})
target_include_directories(compress_mgard PRIVATE include)
target_link_libraries(compress_mgard PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(compress_mgard_sort ${SRC_COMPRESS_MGARD_SORT_FILES})
target_include_directories(compress_mgard_sort PRIVATE include)
target_link_libraries(compress_mgard_sort PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(compress_adios ${SRC_COMPRESS_ADIOS_FILES})
target_include_directories(compress_adios PRIVATE include)
target_link_libraries(compress_adios PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(compress_with_adios_plugin ${SRC_COMPRESS_ADIOS_PLUGIN_FILES})
target_include_directories(compress_with_adios_plugin PRIVATE include)
target_link_libraries(compress_with_adios_plugin PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(decompress_with_adios_plugin ${SRC_DECOMPRESS_ADIOS_PLUGIN_FILES})
target_include_directories(decompress_with_adios_plugin PRIVATE include)
target_link_libraries(decompress_with_adios_plugin PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(evaluate_compression_error ${SRC_EVALUATE_ERROR_FILES})
target_include_directories(evaluate_compression_error PRIVATE include)
target_link_libraries(evaluate_compression_error PRIVATE adios2::cxx11_mpi mgard::mgard MPI::MPI_CXX)

add_executable(compress_zfp ${SRC_COMPRESS_ZFP_FILES})
target_include_directories(compress_zfp PRIVATE include)
target_link_libraries(compress_zfp PRIVATE adios2::cxx11_mpi mgard::mgard zfp::zfp MPI::MPI_CXX)

add_executable(compress_sz ${SRC_COMPRESS_SZ_FILES})
target_include_directories(compress_sz PRIVATE include)
target_link_libraries(compress_sz PRIVATE adios2::cxx11_mpi mgard::mgard zfp::zfp SZ3::SZ3 MPI::MPI_CXX INTERFACE PkgConfig::ZSTD)