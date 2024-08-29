adios2_ROOT="/ccs/home/reshniakv/libs_frontier/adios"
export CMAKE_PREFIX_PATH=$adios2_ROOT:$CMAKE_PREFIX_PATH

mgard_ROOT="/ccs/home/reshniakv/MGARD/install-hip-frontier"
export CMAKE_PREFIX_PATH=$mgard_ROOT:$CMAKE_PREFIX_PATH

cmake -S . -B build
cmake --build build
