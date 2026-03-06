module load PrgEnv-gnu
module load cray-mpich
# module load ums/default
# module load ums002/default
# module load sz3/3.1.7
# module load zstd/1.5.0

adios2_ROOT="/ccs/home/reshniakv/libs_frontier/adios"
export CMAKE_PREFIX_PATH=$adios2_ROOT:$CMAKE_PREFIX_PATH

zfp_ROOT="/ccs/home/reshniakv/libs_frontier/zfp"
export CMAKE_PREFIX_PATH=$zfp_ROOT:$CMAKE_PREFIX_PATH

sz3_ROOT="/ccs/home/reshniakv/libs_frontier/sz3"
export CMAKE_PREFIX_PATH=$sz3_ROOT:$CMAKE_PREFIX_PATH

mgard_ROOT="/ccs/home/reshniakv/MGARD/install-hip-frontier"
export CMAKE_PREFIX_PATH=$mgard_ROOT:$CMAKE_PREFIX_PATH

cmake -S . -B build
cmake --build build

