#!/bin/bash
#PBS -l nodes=2:ppn=28
#PBS -o outnew.log
#PBS -j oe

#cd /home/netweb/vinod_work/breaking_wave/
cd $PBS_O_WORKDIR

export I_MPI_DAPL_DIRECT_COPY_THRESHOLD=655360
export I_MPI_SHM_LMT=shm
export I_MPI_EAGER_THRESHOLD=128000
export MV2_SMP_USE_CMA=0
date
cat $PBS_NODEFILE |uniq > nodes
# qcc -source -grid=octree -D_MPI=1 breaking_iaf_3D_with_column.c.c

mpicc -Wall -std=c99 _breaking_iaf_3D_with_column.c  -o parallel -lm
mpiexec --hostfile nodes -np 56 ./parallel > output.txt

date
