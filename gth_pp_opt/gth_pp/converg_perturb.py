#! /usr/env/bin/python

import os
import subprocess
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op

def run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, \
                        element, method, val_elec_num, python_exe, get_min_index):

  '''
  run_converg_perturb : run convergence perturbation process

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    gth_pp_file : string
      gth_pp_file is the inital guess gth pp.
    cp2k_exe : string
      cp2k_exe is the cp2k exacutable file.
    parallel_exe : string
      parallel_exe is the GNU parallel exacutable file.
    element : string
      element is atom name.
    method : string
      method is self-consistent field method. For DFT, the method is its functional. For HF, the method is HF.
    val_elec_num : int
      val_elec_num is the valence electron number for a atom.
    python_exe : string
      python_exe is the python exacutable file.
    get_min_index : string
      get_min_index is the get_min_index file.
  Returns :
    The restart index where the value is the smallest.
  '''

  #make process_4 directory, and then copy initial guess gth pp file in process_4. 
  #The initial guess gth pp file for process_4 is produced by process_3.
  process_4_dir = ''.join((work_dir, '/process_4'))
  if ( not os.path.exists(process_4_dir) ):
    subprocess.run('mkdir process_4', cwd=work_dir, shell=True)
  subprocess.run('cp %s %s' % (gth_pp_file, process_4_dir), cwd=work_dir, shell=True)

  cmd = "sed -ie '1s/.*/%s GTH-%s-q%d/' GTH-PARAMETER" % (element, method, val_elec_num)
  call.call_simple_shell(process_4_dir, cmd)

  cmd = "ls | grep %s" %('restart')
  restart_num = len(call.call_returns_shell(process_4_dir, cmd))
  if ( restart_num != 0 ):
    cmd = "ls | grep %s" %('bak_1')
    bak_num = len(call.call_returns_shell(process_4_dir, cmd))
    bak_dir = ''.join((process_4_dir, '/bak_', str(bak_num+1)))
    cmd = "mkdir %s" %(''.join(('bak_', str(bak_num+1))))
    call.call_simple_shell(process_4_dir, cmd)
    cmd = "mv restart* %s" %(''.join(('bak_', str(bak_num+1))))
    call.call_simple_shell(process_4_dir, cmd)

  perturb = '''
#! /bin/bash

direc=%s

conv=0.001
weight_standard='2 5 1'
converge_standard='0.003 0.0003 0.003'
weight_perturb_choice_1='30 2 1'
converge_perturb_choice_1='0.0003 0.00003 0.0003'
converge_perturb_choice_2='0.0002 0.00002 0.0002'
converge_perturb_choice_3='0.0001 0.00001 0.0001'
converge_perturb_choice_4='0.03 0.003 0.03'
converge_perturb_choice_5='0.008 0.0008 0.008'
converge_perturb_choice_6='0.007 0.0007 0.007'
converge_perturb_choice_7='0.006 0.0006 0.006'
converge_perturb_choice_8='0.005 0.0005 0.005'
converge_perturb_choice_9='0.004 0.0004 0.004'
converge_perturb_choice_10='0.0003 0.0003 0.0003'
converge_perturb_choice_11='0.003 0.003 0.003'
converge_perturb_choice_12='0.01 0.01 0.001'

m=4
n=4
base_1=4
base_2=1
choice=1
filename=value

$direc/optimize.sh "${weight_standard}" "${converge_perturb_choice_1}" 1
value_1=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" 2`
echo 2 $value_1
echo 2 $value_1 >> $direc/$filename

$direc/optimize.sh "${weight_standard}" "${converge_perturb_choice_1}" 3 new 2
value_2=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" 4 new 3`
echo 4 $value_2
echo 4 $value_2 >> $direc/$filename

for i in {1..2000}
do
if [[ `echo "$(echo "scale=4; $value_1 - $value_2" | bc) > $conv" | bc` == 1 ]]; then
choice=$choice
value_1=$value_2
((j=base_1+(($i-base_2))*2+1))
converge_perturb_new=`eval echo '$'"converge_perturb_choice_$choice"`
$direc/optimize.sh "${weight_standard}" "${converge_perturb_new}" $j old $m
((k=base_1+(($i-base_2))*2+2))
value_2=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" $k old $j`
((m=$m+2))
n=$m
echo $k $value_2
echo $k $value_2 >> $direc/$filename

elif [[ `echo "$(echo "scale=4; $value_1 - $value_2" | bc) < $conv" | bc` == 1 ]]; then
if [ $choice -lt 12 ]; then
((choice=$choice+1))
elif [ $choice == 12 ]; then
choice=1
fi
if [[ `echo "$echo $value_1 < $value_2" | bc` == 1 ]]; then
value_1=$value_1
((j=base_1+(($i-base_2))*2+1))
((n=$n-2))
converge_perturb_new=`eval echo '$'"converge_perturb_choice_$choice"`
$direc/optimize.sh "${weight_standard}" "${converge_perturb_new}" $j new $n
((k=base_1+(($i-base_2))*2+2))
value_2=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" $k new $j`
((n=$n+2))
((m=$m+2))

elif [[ `echo "$echo $value_1 < $value_2" | bc` == 0 ]]; then
value_1=$value_2
((j=base_1+(($i-base_2))*2+1))
converge_perturb_new=`eval echo '$'"converge_perturb_choice_$choice"`
$direc/optimize.sh "${weight_standard}" "${converge_perturb_new}" $j new $m
((k=base_1+(($i-base_2))*2+2))
value_2=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" $k new $j`
((m=$m+2))
n=$m
fi
echo $k $value_2
echo $k $value_2 >> $direc/$filename
fi
done
''' % (process_4_dir)

  optimize = '''
#! /bin/bash

direc=%s
parallel_exe=%s
element=%s
method=%s
elec_num=%d
python_exe=%s
get_index_py=%s

produce() {
x=$1
y=$2
direc=$3
cd $direc/restart$y/step_$x
cp2k_exe=%s
$cp2k_exe atom.inp 1> atom.out 2> atom.err
}
export -f produce

initial_total_step=20
normal_total_step=11
initial_parallel_num=20
normal_parallel_num=11

weight_str=$1
conv_str=$2

weight_arr=(${weight_str// / })
conv_arr=(${conv_str// / })

weight_semi=${weight_arr[0]}
weight_val=${weight_arr[1]}
weight_vir=${weight_arr[2]}

conv_semi=${conv_arr[0]}
conv_val=${conv_arr[1]}
conv_vir=${conv_arr[2]}

i=$3
check_direc=$4
kk=$5

mkdir $direc/restart$i
if [[ $i == 1 || $i == 2 || $check_direc == new ]]; then
for j in $(seq 1 $initial_total_step)
do
mkdir $direc/restart$i/step_$j
cp $direc/../atom.inp $direc/restart$i/step_$j
done
for j in $(seq 1 10)
do
step_size=`echo "scale=6; 0.011-$j*0.001" | bc`
sed -ie '36s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
done
for j in $(seq 11 19)
do
step_size=`echo "scale=6; 0.001-($j-10)*0.0001" | bc`
sed -ie '36s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
done
sed -ie '36s/.*/    STEP_SIZE  0.00009/' $direc/restart$i/step_20/atom.inp
fi
if [[ $i != 1 && $i != 2 && $check_direc == old ]]; then
for j in $(seq 1 $normal_total_step)
do
mkdir $direc/restart$i/step_$j
cp $direc/../atom.inp $direc/restart$i/step_$j
done
fi

sed -ie '45s/.*/     WEIGHT_POT_SEMICORE               '$weight_semi'/' $direc/restart$i/step*/atom.inp
sed -ie '46s/.*/     WEIGHT_POT_VALENCE                '$weight_val'/' $direc/restart$i/step*/atom.inp
sed -ie '47s/.*/     WEIGHT_POT_VIRTUAL                '$weight_vir'/' $direc/restart$i/step*/atom.inp
sed -ie '41s/.*/     TARGET_POT_SEMICORE      [eV]      '$conv_semi'/' $direc/restart$i/step*/atom.inp
sed -ie '42s/.*/     TARGET_POT_VALENCE       [eV]      '$conv_val'/' $direc/restart$i/step*/atom.inp
sed -ie '43s/.*/     TARGET_POT_VIRTUAL       [eV]      '$conv_vir'/' $direc/restart$i/step*/atom.inp
rm $direc/restart$i/step*/*e

if [[ $i == 1 || $kk == 0 ]]; then
for j in $(seq 1 $initial_total_step)
do
sed -ie '29s/.*/    POTENTIAL_FILE_NAME  ..\/..\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done

seq 1 $initial_total_step | $parallel_exe -j $initial_parallel_num produce {} $i $direc

elif [[ $i == 2 || $kk == 1 ]]; then
((k=$i-1))
cd $direc/restart$k
grep -H "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k/kk`
sed -ie '1s/.*/'$element' GTH-'$method'-q'$elec_num'/' step_$m/GTH-PARAMETER
for j in $(seq 1 $initial_total_step)
do
sed -ie '29s/.*/    POTENTIAL_FILE_NAME  ..\/..\/restart'$k'\/step_'$m'\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done

seq 1 $initial_total_step | $parallel_exe -j $initial_parallel_num produce {} $i $direc

elif [[ $i != 1 && $i != 2 ]]; then
k_1=$kk
cd $direc/restart$k_1
grep -H "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k_1/kk`
sed -ie '1s/.*/'$element' GTH-'$method'-q'$elec_num'/' step_$m/GTH-PARAMETER
if [ $check_direc == new ]; then
total_step=$initial_total_step
elif [ $check_direc == old ]; then
total_step=$normal_total_step
fi
for j in $(seq 1 $total_step)
do
sed -ie '29s/.*/    POTENTIAL_FILE_NAME  ..\/..\/restart'$k_1'\/step_'$m'\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done
if [ $check_direc == old ]; then
((k_2=$kk-1))
cd $direc/restart$k_2
grep -H "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k_2/kk`
a=`sed -n '36p' step_$m/atom.inp`
step_size_pre=`echo $a | tr -cd "[0-9,.]"`
for j in $(seq 1 $total_step)
do
scalling=`echo "scale=6; 0.4+0.1*$j" | bc`
step_size=`echo "scale=6; $step_size_pre*$scalling" | bc`
sed -ie '36s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
done
fi

if [ $check_direc == new ]; then
parallel_num=$initial_parallel_num
elif [ $check_direc == old ]; then
parallel_num=$normal_parallel_num
fi
seq 1 $total_step | $parallel_exe -j $parallel_num produce {} $i $direc
fi

cd $direc/restart$i
for z in {1..20}; do grep -H "Final value of function" step_$z/atom.out; done > kk
m=`$python_exe $get_index_py $direc/restart$i/kk`
str=`grep -H "Final value of function" step_$m/atom.out`
value=`echo ${str:20:100} | tr -cd "[0-9,.]"`
echo $value
cd $direc
''' %(process_4_dir, parallel_exe, element, method, val_elec_num, python_exe, get_min_index, cp2k_exe)

  #run shell script, and then get returns
  optimize_file = ''.join((process_4_dir, '/optimize.sh'))
  with open(optimize_file, 'w') as f:
    f.write(optimize)

  perturb_file = ''.join((process_4_dir, '/perturb.sh'))
  with open(perturb_file, 'w') as f:
    f.write(perturb)

  subprocess.run('chmod +x perturb.sh', cwd=process_4_dir, shell=True)
  subprocess.run('chmod +x optimize.sh', cwd=process_4_dir, shell=True)
  p = subprocess.Popen("bash -c './perturb.sh'", cwd=process_4_dir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  #kill the shell process if the result does not decay.
  restart_index = []
  value = []
  for line in iter(p.stdout.readline, 'b'):
    line = line.rstrip().decode('utf8')
    if ( line != '' ):
      line_split = data_op.split_str(line, ' ')
      if ( len(line_split) == 2 ):
        restart_index.append(int(line_split[0]))
        value.append(float(line_split[1]))
    if ( subprocess.Popen.poll(p) is not None ):
      if ( line == '' ):
        break
    if ( len(value) > 40 ):
      min_value = min(value)
      if ( len(line_split) == 2 ):
        if ( min_value < float(line_split[1]) ):
          call.kills(p.pid)

  min_value = min(value)
  min_value_index = value.index(min_value)

  return restart_index[min_value_index]

if __name__ == '__main__':

  from CP2K_kit.gth_pp_opt.gth_pp import converg_perturb

  work_dir = '/home/lujunbo/code/github/CP2K_kit/example/gth_pp_opt'
  gth_pp_file = '/home/lujunbo/code/github/CP2K_kit/example/gth_pp_opt/process_3/restart2/step_12/GTH-PARAMETER'
  cp2k_exe = '/home/lujunbo/bin/cp2k/cp2k-6.1-sopt/cp2k-6.1-Linux-x86_64.sopt'
  parallel_exe = '/home/lujunbo/bin/parallel/bin/parallel'
  element = 'Fm'
  method = 'PBE'
  val_elec_num = 22
  python_exe = 'python3.8'
  get_min_index = '/home/lujunbo/code/github/CP2K_kit/gth_pp_opt/gth_pp/get_index.py'
  converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, \
                                      element, method, val_elec_num, python_exe, get_min_index)
