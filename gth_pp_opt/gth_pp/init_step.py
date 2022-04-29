#! /usr/env/bin python

import subprocess
from CP2K_kit.tools import file_tools

def gen_init_step(work_dir, gth_pp_file):

  '''
  gen_init_step : generate initial gth pp optimization directory.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    gth_pp_file : string
      gth_pp_file is the inital guess gth pp.
  Returns :
    none
  '''

  atom_inp_file = ''.join((work_dir, '/atom.inp'))
  line_num = file_tools.grep_line_num('STEP_SIZE', atom_inp_file, work_dir)

  gen_direc = '''
#! /bin/bash

direc=%s
line_num=%d
for i in {1..129}
do
if [ ! -f $direc/process_1/step_$i ]; then
mkdir $direc/process_1/step_$i
fi
cp $direc/atom.inp $direc/process_1/step_$i
done

for i in {1..100}
do
step_size=`echo "scale=6; 10.1-$i*0.1" | bc`
sed -ie ''${line_num}'s/.*/    STEP_SIZE  '$step_size'/' $direc/process_1/step_$i/atom.inp
done

for i in {101..110}
do
step_size=`echo "scale=6; 0.1-($i-100)*0.01" | bc`
sed -ie ''${line_num}'s/.*/    STEP_SIZE  '$step_size'/' $direc/process_1/step_$i/atom.inp
done

for i in {111..120}
do
step_size=`echo "scale=6; 0.01-($i-110)*0.001" | bc`
sed -ie ''${line_num}'s/.*/    STEP_SIZE  '$step_size'/' $direc/process_1/step_$i/atom.inp
done

for i in {121..129}
do
step_size=`echo "scale=6; 0.001-($i-120)*0.0001" | bc`
sed -ie ''${line_num}'s/.*/    STEP_SIZE  '$step_size'/' $direc/process_1/step_$i/atom.inp
done
''' % (work_dir, line_num[0])

  #make process_1 directory, and copy initial gth file in process_1 directory
  #The initial gth file in process_1 is defined by user.
  process_1_dir = ''.join((work_dir, '/process_1'))
  subprocess.run("mkdir %s" % (process_1_dir), cwd=work_dir, shell=True)
  process_1_dir_pp = ''.join((process_1_dir, '/GTH-PARAMETER'))
  subprocess.run("cp %s %s" % (gth_pp_file, process_1_dir_pp), cwd=work_dir, shell=True)

  init_step_file = ''.join((work_dir, '/process_1/init_step.sh'))
  with open(init_step_file, 'w') as f:
    f.write(gen_direc)

  process_1_dir = ''.join((work_dir, '/process_1'))
  subprocess.run('chmod +x init_step.sh', cwd=process_1_dir, shell=True)
  subprocess.run("bash -c './init_step.sh'", cwd=process_1_dir, shell=True)

def run_init_step(work_dir, cp2k_exe, parallel_exe, element, val_elec_num, method, parallel_num, start, end):

  '''
  run_init_step : perform gth pp optimization for initial steps. There are 129 steps.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    cp2k_exe : string
      cp2k_exe is the cp2k exacutable file.
    parallel_exe : string
      parallel_exe is the GNU parallel exacutable file.
    element : string
      element is atom name.
    val_elec_num : int
      val_elec_num is the valence electron number for a atom.
    method : string
      method is self-consistent field method. For DFT, the method is its functional. For HF, the method is HF.
    parallel_num : int
      parallel_num is the number of processes for parallel program.
    start : int
      start is the starting step.
    end : int
      end is the ending step.
  Returns:
    none
  '''

  run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

produce() {
x=$1
direc=$2
element=%s
elec_num=%d
method=%s
cp2k_exe=%s
cd $direc/process_1/step_$x
$cp2k_exe atom.inp 1> atom.out 2> atom.err
cd $direc/process_1
}

export -f produce

seq $run_start $run_end | $parallel_exe -j $parallel_num produce {} $direc
''' %(work_dir, parallel_num, start, end, parallel_exe, element, val_elec_num, method, cp2k_exe)

  process_1_dir = ''.join((work_dir, '/process_1'))
  run_file = ''.join((work_dir, '/process_1/run.sh'))
  with open(run_file, 'w') as f:
    f.write(run)

  subprocess.run('chmod +x run.sh', cwd=process_1_dir, shell=True)
  subprocess.run("bash -c './run.sh'", cwd=process_1_dir, shell=True)
