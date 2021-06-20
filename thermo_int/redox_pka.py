#! /usr/env/bin python

import subprocess

def slow_growth(inp_file, start, end, restart, proc_num, cp2k_exe, env_cp2k):

slow_growth = '''
#! /bin/sh

project=%s
initial_file_name=$s
initial_step=%d
restart_step=$1
end_step=%d
scope=1
process_num=%d
cp2k_exe=%s

increment=`echo "scale=6; $scope/($end_step-$initial_step)" | bc`

for i in $(seq $restart_step $end_step)
do
if [ $i == $restart_step ];then
cp restart_bak_file $project'-1.restart'
fi
if [ $i -gt $initial_step ];then
((j=$i-$initial_step))
a=`echo "scale=6; $increment*$j" | bc`
line=`grep -n "VALUE" $project'-1.restart'`
n_1=${line:0:3}
n_2=${line:0:2}
expr $n_1 "+" 10 &> /dev/null
if [ $? -eq 0 ];then
n=$n_1
else
n=$n_2
fi
if  [ ! -f  $project'-1.restart'  ]; then
((m=$i-1))
echo $m
exit
else
sed -ie ''$n's/.*/       VALUES     '$a'/' $project'-1.restart'
mv $project'-1.restart' input.inp
rm $project'-1.restarte'
fi
mpirun -np $process_num $cp2k_exe input.inp > cp2k.out 2> cp2k.err
rm input.inp
cp $project'-1.restart' restart_bak_file
else
mpirun -np $process_num $cp2k_exe $initial_file_name > cp2k.out 2> cp2k.err
fi
done
''' % (inp_file, start, end, proc_num, cp2k_exe)

  slow_growth_file = ''.join((work_dir, '/slow_growth.sh'))
  with open(slow_growth_file, 'w') as f:
    f.write(slow_growth)

  subprocess.run('chmod +x slow_growth.sh', cwd=work_dir, shell=True)

run = '''
#! /bin/bash

a=`./slow_growth_ti.sh %d`
for i in {1..100}
do
if [ ! -n "$a" ]; then
exit
else
a=`./slow_growth_ti.sh $a`
fi
done
''' % (restart)

  run_file = ''.join((work_dir, '/run.sh'))
  with open(run_file, 'w') as f:
    f.write(run)

  subprocess.run('chmod +x run.sh', cwd=work_dir, shell=True)
  subprocess.run("bash -c './run.sh'", cwd=work_dir, shell=True, env=env_cp2k)

def run_slow_growth(work_dir, input_file):

  if ( 'inp_file' in redox_pka_dic.keys() ):
    inp_file = redox_pka_dic['inp_file']
  else:
    print ('Cannot find input file, please set inp_file')
    exit()

  if ( 'start' in redox_pka_dic.keys() ):
    start = int(redox_pka_dic['start'])
  else:
    print ('Cannot find starting step, please set start')
    exit()

  if ( 'end' in redox_pka_dic.keys() ):
    end = int(redox_pka_dic['end'])
  else:
    print ('Cannot find endding step, please set end')
    exit()

  if ( 'restart' in redox_pka_dic.keys() ):
    restart = int(redox_pka_dic['restart'])
  else:
    print ('Cannot find restarting step, please set restart')
    exit()

  if ( 'proc_num' in redox_pka_dic.keys() ):
    proc_num = int(redox_pka_dic['proc_num'])
  else:
    print ('Cannot find process number, please set proc_num')
    exit()

  if ( 'cp2k_exe' in redox_pka_dic.keys() ):
    cp2k_exe = redox_pka_dic['cp2k_exe']
  else:
    print ('Cannot find cp2k executable file, please set cp2k_exe')
    exit()

  if ( 'mpi_dir' in redox_pka_dic.keys() )
    mpi_dir = redox_pka_dic['mpi_dir']
  else:
    print ('Cannot find mpi directory, please set mpi_dir')
    exit()

  env_cp2k = {}
  env_cp2k["PATH"] = ''.join((mpi_dir, '/bin'))
  env_cp2k["LD_LIBRARY_PATH"] = ''.join((mpi_dir, '/lib'))

  slow_growth(inp_file, start, end, restart, proc_num, cp2k_exe, env_cp2k)



