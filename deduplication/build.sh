module load PrgEnv-gnu
module load cray-mpich

adios2_ROOT="/ccs/home/reshniakv/libs_frontier/adios"
export CMAKE_PREFIX_PATH=$adios2_ROOT:$CMAKE_PREFIX_PATH

cmake -S . -B build
cmake --build build

