#!/bin/bash
#PBS -q long
#PBS -l nodes=1:ppn=24
#PBS -j oe
#PBS -V
#PBS -N UO22+_q24
#PBS -o output.o
#PBS -e error.e
NP=`cat $PBS_NODEFILE | wc -l`
cd $PBS_O_WORKDIR
##### definition of job name ################
JOB_NAME=cal
#echo -n "start time  " > time

echo "job ${JOB_NAME} starts at `date`" >${JOB_NAME}.out
echo "running on the following nodes, with $NP processes in total" >>${JOB_NAME}.out
cat $PBS_NODEFILE | sort | uniq -c >>${JOB_NAME}.out

#--------------------- cp2k calculation-------------------------------
for file_a in ${PBS_O_WORKDIR}/*.inp; do
echo "Start Time:" `date` > time
date_start=`date "+%s"`
FDNAME=`basename $PBS_O_WORKDIR`
ulimit -s unlimited
export PATH=/home/lujunbo/bin/cp2k-8.2/exe/local:$PATH
source /home/lujunbo/bin/cp2k-8.2/tools/toolchain/install/setup
mpirun -np 24 cp2k.popt cp2k-aimd.inp 1>> qmmm.out 2>> qmmm.err
#############################
echo "END Time:" `date` >> time
date_end=`date "+%s"`
hour=`awk -v y=$date_start -v x=$date_end 'BEGIN {printf "%.2f\n",(x-y)/3600.0}'`
echo "Runing Time(h):" $hour "(h)">> time
done
#---------------------------------------------------------------------
