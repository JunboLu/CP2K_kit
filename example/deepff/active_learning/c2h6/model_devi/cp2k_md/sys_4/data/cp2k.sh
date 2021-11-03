#!/bin/bash
#PBS -q long
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -V
#PBS -N deepmd
#PBS -o output.o
#PBS -e error.e
cd $PBS_O_WORKDIR

source ~/profile/cp2k_8.2.sh

for i in {2748..3000}
do
cd /home/lujunbo/WORK/deepmd/c2h6/mtd/train_4_sys/cp2k_md/sys_4/data/task_$i 
mpirun -np 24 cp2k.popt input.inp 1> qmmm.out 2> qmmm.err
done
