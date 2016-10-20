cmake_minimum_required(VERSION 3.5)
project(fast_protein_cluster)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    test/1af7_.url
    test/cluster_output.agglomeration.history
    test/cluster_output.cluster.stats
    test/cluster_output.clusters
    test/cluster_output.density
    test/list
    test/list.10
    test/list.100
    test/list.1000
    test/list.10000
    test/list.30
    test/list.300
    test/output.agglomeration.history
    test/output.cluster.stats
    test/output.clusters
    test/README
    test/test_rmsd.sh
    test/test_rmsd_avx.sh
    test/test_rmsd_density.sh
    test/test_rmsd_density_avx.sh
    test/test_rmsd_density_gpu.sh
    test/test_rmsd_density_sse2.sh
    test/test_rmsd_density_sse3.sh
    test/test_rmsd_gpu.sh
    test/test_rmsd_sse2.sh
    test/test_rmsd_sse3.sh
    test/test_rmsd_tmscore.sh
    test/test_rmsd_tmscore_avx.sh
    test/test_rmsd_tmscore_gpu.sh
    test/test_rmsd_tmscore_sse2.sh
    test/test_rmsd_tmscore_sse3.sh
    test/test_tmscore.sh
    test/test_tmscore_avx.sh
    test/test_tmscore_gpu.sh
    test/test_tmscore_sse2.sh
    test/test_tmscore_sse3.sh
    CHANGE_LOG
    data_structures.h
    defines.h
    error_handlers.c
    error_handlers.h
    fast_protein_cluster.cpp
    libcluster.hpp
    libtmscore_cpu.cpp
    libtmscore_gpu.cpp
    LICENSE.txt
    lite.h
    Makefile
    my_sort.h
    OPTIONS.txt
    rand31pmc.cpp
    rand31pmc.h
    README
    README_GPU
    rmsd.cl
    tmscore.cl
    torsions.h)

add_executable(fast_protein_cluster ${SOURCE_FILES})