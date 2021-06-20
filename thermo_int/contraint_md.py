#! /usr/env/bin python

def constraint_md():

incre_point = '''
#! /bin/bash

direc=%s
project=%s
start_point=%f
increment=%f
initial_step=%d
restart_step=%d
end_step=%d
total_md_step=%d
process_num=%d
cp2k_exe=%f

inp_file=%s

for i in $(seq $restart_step $end_step)
do
mkdir $direc/point_$i

if [ $i == $initial_step ]; then
cp $inp_file $project-RESTART.wfn $direc/point_$i
cd $direc/point_$i
line=`grep -n "TARGET" input.inp`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
sed -ie ''$n's/.*/      TARGET [angstrom] '$start_point'/' input.inp
mpirun -np $process_num $cp2k_direc_exe input.inp > qmmm.out 2> qmmm.err
cd ..
fi

if [[ $i != $initial_step && $i == $restart_step ]]; then
cd $direc/point_$i
a_1=`grep "STEP_START_VAL" $project-1.restart`
a_2=`grep "STEP_START_VAL" input.inp`
b_1=`echo $a_1 | tr -cd "[0-9]"`
b_2=`echo $a_2 | tr -cd "[0-9]"`
((num=b_1-b_2))
((remain_num=total_md_step-num))
line=`grep -n " STEPS " $project-1.restart`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
sed -ie ''$n's/.*/     STEPS   '$remain_num'/' $project-1.restart
rm $project-1.restarte
cp $project-RESTART.wfn.bak-1 $project-RESTART.wfn.bak

a=`ls | grep "qmmm" | grep "out" | grep "bak"`
b=`echo $a | tr -cd "[0-9]"`
if [ -z $b ]; then
mpirun -np $process_num $cp2k_exe $project-1.restart > qmmm_bak_1.out 2> qmmm_bak_1.err
else
((c=b+1))
mpirun -np $process_num $cp2k_exe $project-1.restart > qmmm_bak_$c.out 2> qmmm_bak_$c.err
fi
cd ..
fi

if [[ $i != $initial_step && $i != $restart_step ]]; then
((k=i-1))
cd $direc/point_$k
cp $project-1.restart $project-RESTART.wfn $direc/point_$i
cd $direc/point_$i
a=`echo "scale=6; $start_point+$increment*$k" | bc`

line=`grep -n "TARGET" $project-1.restart`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
sed -ie ''$n's/.*/      TARGET [angstrom] '$a'/' $project-1.restart

line=`grep -n " STEPS " $project-1.restart`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
sed -ie ''$n's/.*/     STEPS   '$total_md_step'/' $project-1.restart

rm $project-1.restarte*
mv $project-1.restart input.inp
mpirun -np $process_num $cp2k_exe input.inp > qmmm.out 2> qmmm.err
fi

done
'''

slow_growth = '''
#! /bin/bash

mpirun -np 96 /home/lujunbo/bin/code_block/cp2k-6.1.0/arch/Linux-x86-64-intel/cp2k.popt UO22+_aimd-1.restart > qmmm.out 2> qmmm.err

start_point=2.9140959719616807

for i in {1..15131}
do
a=`echo "scale=6; $start_point+0.0001*$i" | bc`
line=`grep -n "TARGET" UO22+_aimd-1.restart`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
sed -ie ''$n's/.*/      TARGET [angstrom] '$a'/' UO22+_aimd-1.restart
rm UO22+_aimd-1.restarte
mpirun -np 96 /home/lujunbo/bin/code_block/cp2k-6.1.0/arch/Linux-x86-64-intel/cp2k.popt UO22+_aimd-1.restart > qmmm.out 2> qmmm.err
done
'''
