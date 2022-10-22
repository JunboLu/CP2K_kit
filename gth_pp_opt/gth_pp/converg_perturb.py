#! /usr/env/bin/python

import os
import time
import tempfile
import subprocess
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import file_tools

def run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                        method, val_elec_num, python_exe, get_min_index, weight_2,
                        converge_perturb_choice_1, converge_perturb_choice_2, \
                        converge_perturb_choice_3, converge_perturb_choice_4, \
                        converge_perturb_choice_5, converge_perturb_choice_6, \
                        converge_perturb_choice_7, converge_perturb_choice_8, \
                        converge_perturb_choice_9, converge_perturb_choice_10, \
                        converge_perturb_choice_11, converge_perturb_choice_12, \
                        target_semi, target_val, target_vir, consider_charge, max_cycle=2000):

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
    weight_2: 1-d float list
      weight_2 is the initial weight 2.
  Returns :
    The restart index where the value is the smallest.
  '''

  #make process_4 directory, and then copy initial guess gth pp file in process_4. 
  #The initial guess gth pp file for process_4 is produced by process_3.
  if ( max_cycle != 2000 ):
    process_4_dir = ''.join((work_dir, '/process_mix'))
  else:
    process_4_dir = ''.join((work_dir, '/process_4'))
  if ( (not os.path.exists(process_4_dir)) and max_cycle != 2000 ):
    subprocess.run('mkdir process_mix', cwd=work_dir, shell=True)
  if ( (not os.path.exists(process_4_dir)) and max_cycle == 2000 ):
    subprocess.run('mkdir process_4', cwd=work_dir, shell=True)
  subprocess.run('cp %s %s' % (gth_pp_file, process_4_dir), cwd=work_dir, shell=True)

  cmd = "sed -ie '1s/.*/%s GTH-%s-q%d/' GTH-PARAMETER" % (element, method, val_elec_num)
  call.call_simple_shell(process_4_dir, cmd)

  cmd = "ls | grep %s" %('restart')
  restart_num = len(call.call_returns_shell(process_4_dir, cmd))
  if ( restart_num != 0 ):
    cmd = "ls | grep %s" %('bak_')
    bak_num = len(call.call_returns_shell(process_4_dir, cmd))
    bak_dir = ''.join((process_4_dir, '/bak_', str(bak_num+1)))
    cmd = "mkdir %s" %(''.join(('bak_', str(bak_num+1))))
    call.call_simple_shell(process_4_dir, cmd)
    cmd = "mv restart* optimize.sh perturb.sh %s" %(''.join(('bak_', str(bak_num+1))))
    call.call_simple_shell(process_4_dir, cmd)

  converge_perturb_choice_1_str = data_op.comb_list_2_str(converge_perturb_choice_1, ' ')
  converge_perturb_choice_2_str = data_op.comb_list_2_str(converge_perturb_choice_2, ' ')
  converge_perturb_choice_3_str = data_op.comb_list_2_str(converge_perturb_choice_3, ' ')
  converge_perturb_choice_4_str = data_op.comb_list_2_str(converge_perturb_choice_4, ' ')
  converge_perturb_choice_5_str = data_op.comb_list_2_str(converge_perturb_choice_5, ' ')
  converge_perturb_choice_6_str = data_op.comb_list_2_str(converge_perturb_choice_6, ' ')
  converge_perturb_choice_7_str = data_op.comb_list_2_str(converge_perturb_choice_7, ' ')
  converge_perturb_choice_8_str = data_op.comb_list_2_str(converge_perturb_choice_8, ' ')
  converge_perturb_choice_9_str = data_op.comb_list_2_str(converge_perturb_choice_9, ' ')
  converge_perturb_choice_10_str = data_op.comb_list_2_str(converge_perturb_choice_10, ' ')
  converge_perturb_choice_11_str = data_op.comb_list_2_str(converge_perturb_choice_11, ' ')
  converge_perturb_choice_12_str = data_op.comb_list_2_str(converge_perturb_choice_12, ' ')

  atom_inp_file = ''.join((work_dir, '/atom.inp'))
  line_step_size = file_tools.grep_line_num('STEP_SIZE', atom_inp_file, work_dir)
  line_wt_semi = file_tools.grep_line_num('WEIGHT_POT_SEMICORE', atom_inp_file, work_dir)
  line_wt_val = file_tools.grep_line_num('WEIGHT_POT_VALENCE', atom_inp_file, work_dir)
  line_wt_vir = file_tools.grep_line_num('WEIGHT_POT_VIRTUAL', atom_inp_file, work_dir)
  line_tg_semi = file_tools.grep_line_num('TARGET_POT_SEMICORE', atom_inp_file, work_dir)
  line_tg_val = file_tools.grep_line_num('TARGET_POT_VALENCE', atom_inp_file, work_dir)
  line_tg_vir = file_tools.grep_line_num('TARGET_POT_VIRTUAL', atom_inp_file, work_dir)
  line_pot_file = file_tools.grep_line_num('POTENTIAL_FILE_NAME', atom_inp_file, work_dir)

  perturb = '''
#! /bin/bash

direc=%s

conv=0.001
weight_standard='%f %f %f'
converge_standard='%f %f %f'
weight_perturb_choice_1='30 2 1'
converge_perturb_choice_1='%s'
converge_perturb_choice_2='%s'
converge_perturb_choice_3='%s'
converge_perturb_choice_4='%s'
converge_perturb_choice_5='%s'
converge_perturb_choice_6='%s'
converge_perturb_choice_7='%s'
converge_perturb_choice_8='%s'
converge_perturb_choice_9='%s'
converge_perturb_choice_10='%s'
converge_perturb_choice_11='%s'
converge_perturb_choice_12='%s'

m=4
n=4
base_1=4
base_2=1
choice=1
filename=value

$direc/optimize.sh "${weight_standard}" "${converge_perturb_choice_1}" 1
value_1=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" 2`
echo 2 $value_1 $choice
echo 2 $value_1 >> $direc/$filename

$direc/optimize.sh "${weight_standard}" "${converge_perturb_choice_1}" 3 new 2
value_2=`$direc/optimize.sh "${weight_standard}" "${converge_standard}" 4 new 3`
echo 4 $value_2 $choice
echo 4 $value_2 >> $direc/$filename

for i in {1..%d}
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
echo $k $value_2 $choice
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
echo $k $value_2 $choice
echo $k $value_2 >> $direc/$filename
fi
done
''' % (process_4_dir, weight_2[0], weight_2[1], weight_2[2], target_semi, target_val, target_vir, \
       converge_perturb_choice_1_str, converge_perturb_choice_2_str, converge_perturb_choice_3_str, \
       converge_perturb_choice_4_str, converge_perturb_choice_5_str, converge_perturb_choice_6_str, \
       converge_perturb_choice_7_str, converge_perturb_choice_8_str, converge_perturb_choice_9_str, \
       converge_perturb_choice_10_str, converge_perturb_choice_11_str, converge_perturb_choice_12_str, max_cycle)

  optimize = '''
#! /bin/bash

direc=%s
parallel_exe=%s
element=%s
method=%s
elec_num=%d
python_exe=%s
get_index_py=%s
consider_charge=%s

line_step_size=%d
line_wt_semi=%d
line_wt_val=%d
line_wt_vir=%d
line_tg_semi=%d
line_tg_val=%d
line_tg_vir=%d
line_pot_file=%d

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
sed -ie ''${line_step_size}'s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
done
for j in $(seq 11 19)
do
step_size=`echo "scale=6; 0.001-($j-10)*0.0001" | bc`
sed -ie ''${line_step_size}'s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
done
sed -ie ''${line_step_size}'s/.*/    STEP_SIZE  0.00009/' $direc/restart$i/step_20/atom.inp
fi
if [[ $i != 1 && $i != 2 && $check_direc == old ]]; then
for j in $(seq 1 $normal_total_step)
do
mkdir $direc/restart$i/step_$j
cp $direc/../atom.inp $direc/restart$i/step_$j
done
fi

sed -ie ''${line_wt_semi}'s/.*/     WEIGHT_POT_SEMICORE               '$weight_semi'/' $direc/restart$i/step*/atom.inp
sed -ie ''${line_wt_val}'s/.*/     WEIGHT_POT_VALENCE                '$weight_val'/' $direc/restart$i/step*/atom.inp
sed -ie ''${line_wt_vir}'s/.*/     WEIGHT_POT_VIRTUAL                '$weight_vir'/' $direc/restart$i/step*/atom.inp
sed -ie ''${line_tg_semi}'s/.*/     TARGET_POT_SEMICORE      [eV]      '$conv_semi'/' $direc/restart$i/step*/atom.inp
sed -ie ''${line_tg_val}'s/.*/     TARGET_POT_VALENCE       [eV]      '$conv_val'/' $direc/restart$i/step*/atom.inp
sed -ie ''${line_tg_vir}'s/.*/     TARGET_POT_VIRTUAL       [eV]      '$conv_vir'/' $direc/restart$i/step*/atom.inp
rm $direc/restart$i/step*/*e

if [[ $i == 1 || $kk == 0 ]]; then
for j in $(seq 1 $initial_total_step)
do
sed -ie ''${line_pot_file}'s/.*/    POTENTIAL_FILE_NAME  ..\/..\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done

seq 1 $initial_total_step | $parallel_exe -j $initial_parallel_num produce {} $i $direc

elif [[ $i == 2 || $kk == 1 ]]; then
((k=$i-1))
cd $direc/restart$k
grep -H --text "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k $consider_charge`
sed -ie '1s/.*/'$element' GTH-'$method'-q'$elec_num'/' step_$m/GTH-PARAMETER
for j in $(seq 1 $initial_total_step)
do
sed -ie ''${line_pot_file}'s/.*/    POTENTIAL_FILE_NAME  ..\/..\/restart'$k'\/step_'$m'\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done

seq 1 $initial_total_step | $parallel_exe -j $initial_parallel_num produce {} $i $direc

elif [[ $i != 1 && $i != 2 ]]; then
k_1=$kk
cd $direc/restart$k_1
grep -H --text "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k_1 $consider_charge`
sed -ie '1s/.*/'$element' GTH-'$method'-q'$elec_num'/' step_$m/GTH-PARAMETER
if [ $check_direc == new ]; then
total_step=$initial_total_step
elif [ $check_direc == old ]; then
total_step=$normal_total_step
fi
for j in $(seq 1 $total_step)
do
sed -ie ''${line_pot_file}'s/.*/    POTENTIAL_FILE_NAME  ..\/..\/restart'$k_1'\/step_'$m'\/GTH-PARAMETER/' $direc/restart$i/step_$j/atom.inp
done
if [ $check_direc == old ]; then
((k_2=$kk-1))
cd $direc/restart$k_2
grep -H --text "Final value of function" step*/atom.out > kk
m=`$python_exe $get_index_py $direc/restart$k_2 $consider_charge`
a=`sed -n ''${line_step_size}'p' step_$m/atom.inp`
step_size_pre=`echo $a | tr -cd "[0-9,.]"`
for j in $(seq 1 $total_step)
do
scalling=`echo "scale=6; 0.4+0.1*$j" | bc`
step_size=`echo "scale=6; $step_size_pre*$scalling" | bc`
sed -ie ''${line_step_size}'s/.*/    STEP_SIZE  '$step_size'/' $direc/restart$i/step_$j/atom.inp
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
for z in {1..20}; do grep -H --text "Final value of function" step_$z/atom.out; done > kk
m=`$python_exe $get_index_py $direc/restart$i $consider_charge`
str=`grep -H --text "Final value of function" step_$m/atom.out`
value=`echo ${str:20:100} | tr -cd "[0-9,.]"`
echo $value
cd $direc
''' %(process_4_dir, parallel_exe, element, method, val_elec_num, python_exe, get_min_index, \
      consider_charge, line_step_size[0], line_wt_semi[0], line_wt_val[0], line_wt_vir[0], \
      line_tg_semi[0], line_tg_val[0], line_tg_vir[0], line_pot_file[0],cp2k_exe)

  #run shell script, and then get returns
  optimize_file = ''.join((process_4_dir, '/optimize.sh'))
  with open(optimize_file, 'w') as f:
    f.write(optimize)

  perturb_file = ''.join((process_4_dir, '/perturb.sh'))
  with open(perturb_file, 'w') as f:
    f.write(perturb)

  subprocess.run('chmod +x perturb.sh', cwd=process_4_dir, shell=True)
  subprocess.run('chmod +x optimize.sh', cwd=process_4_dir, shell=True)
  out_temp=tempfile.TemporaryFile(mode='w+')
  fileno=out_temp.fileno()
  p = subprocess.Popen("bash -c './perturb.sh'", cwd=process_4_dir, shell=True, stdout=fileno, stderr=fileno)

  #kill the shell process if the result does not decay.
  while p.poll() is None:
    time.sleep(30)
    out_temp.seek(0)
    rt = out_temp.read()
    rt_list = rt.strip().split('\n')
    value = []
    restart_index = []
    choice = []
    for line in iter(rt_list):
      if ( line != '' ):
        line_split = data_op.split_str(line, ' ')
        if ( len(line_split) == 3 and data_op.eval_str(line_split[0]) == 1 and \
             data_op.eval_str(line_split[1]) == 2 and \
             data_op.eval_str(line_split[2]) == 1 ):
          restart_index.append(int(line_split[0]))
          value.append(float(line_split[1]))
          choice.append(int(line_split[2]))
    if ( 12 in choice ):
      index_11 = []
      index_12 = []
      for index, choice_value in enumerate(choice):
        if ( choice_value == 11 ):
          index_11.append(index)
        elif ( choice_value == 12 ):
          index_12.append(index)
      if ( index_11[len(index_11)-1] > index_12[len(index_12)-1] ):
        min_value = min(value[0:(len(value)-1)])
        final_value = value[len(value)-1]
        if ( min_value < final_value or abs(min_value-final_value) < 0.0001 ):
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
