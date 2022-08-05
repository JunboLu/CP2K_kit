#! /usr/env/bin python

import os
import math
import subprocess
import numpy as np
from CP2K_kit.tools import *
from CP2K_kit.deepff import process
from CP2K_kit.deepff import gen_lammps_task

def lmpmd_single(lmp_dir, sys_index, task_index, lmp_exe, lmp_path, mpi_path, \
                 lmp_mpi_num_per_job, lmp_omp_num_per_job, device):

  '''
  lmpmd_single: run single lammps molecular dynamics calculation.

  Args:
    lmp_dir: string
      lmp_dir is the directory of lammps calculation.
    sys_index: int
      sys_index is the index of system.
    task_index: int
      task_index is the index of task.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    lmp_path: string
      lmp_path is the path of lammps.
    mpi_path: string
      mpi_path is the path of mpi.
    lmp_mpi_num_per_job: int
      lmp_mpi_num_per_job is the mpi number for each lammps job.
    lmp_omp_num_per_job: int
      lmp_omp_num_per_job is the openmp number for each lammps job.
    device: 2-d int list
      device is the name of gpu devices for all nodes.
  Returns:
    none
  '''

  lmp_sys_task_dir = ''.join((lmp_dir, '/sys_', str(sys_index), '/task_', str(task_index)))
  file_label = 0
  while True:
    log_file_name = ''.join(('lammps', str(file_label), '.out'))
    log_file_name_abs = ''.join((lmp_sys_task_dir, '/', log_file_name))
    dump_file_name = ''.join(('atom', str(file_label), '.dump'))
    dump_file_name_abs = ''.join((lmp_sys_task_dir, '/', dump_file_name))
    if ( os.path.exists(log_file_name_abs) and os.path.exists(dump_file_name_abs) and os.path.getsize(dump_file_name_abs) != 0 ):
      file_label = file_label+1
    else:
      break

  if ( len(device[0]) == 0 ):
    run = '''
#! /bin/bash

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=%d

mpirun -np %d %s < ./md_in.lammps 1> %s 2> lammps.err
''' %(lmp_path, mpi_path, lmp_omp_num_per_job, lmp_mpi_num_per_job, lmp_exe, log_file_name)

  else:
    device_str=data_op.comb_list_2_str(device[0], ',')
    run = '''
#! /bin/bash

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

export CUDA_VISIBLE_DEVICES=%s
export OMP_NUM_THREADS=%d

mpirun -np %d %s < ./md_in.lammps 1> %s 2> lammps.err
''' %(lmp_path, mpi_path, device_str, lmp_omp_num_per_job, lmp_mpi_num_per_job, lmp_exe, log_file_name)

  run_file_name_abs = ''.join((lmp_sys_task_dir, '/run.sh'))
  with open(run_file_name_abs, 'w') as f:
    f.write(run)

  subprocess.run('chmod +x run.sh', cwd=lmp_sys_task_dir, shell=True)
  try:
    subprocess.run("bash -c './run.sh'", cwd=lmp_sys_task_dir, shell=True)
  except subprocess.CalledProcessError as err:
    log_info.log_error('Running error: %s command running error in %s' %(err.cmd, lmp_sys_task_dir))
    exit()

  gen_lammps_task.combine_frag_traj_file(lmp_sys_task_dir)

def lmpmd_parallel(lmp_dir, lmp_path, mpi_path, lmp_exe, parallel_exe, sys_index_str, \
                   task_index_str, mpi_num_str, device_num_str, device_id_start_str, \
                   lmp_omp_num_per_job, lmp_job_per_node, proc_num_per_node, ssh, host):

  '''
  lmpmd_parallel: run lammps molecular dynamics calculation in parallel.

  Args:
    lmp_dir: string
      lmp_dir is the directory of lammps calculation.
    lmp_path: string
      lmp_path is the path of lammps.
    mpi_path: string
      mpi_path is the path of mpi.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    parallel_exe: string
      parallel_exe is the parallel executable file.
    sys_index_str: string
      sys_index_str is the string containing systems index.
    task_index_str: string
      task_index_str is the string containing tasks index.
    mpi_num_str: string
      mpi_num_str is the string containing mpi number.
    device_num_str: string
      device_num_str is the string containing gpu devices number.
    device_id_start_str: string
      device_id_start_str is the string containing staring gpu device id.
    lmp_omp_num_per_job: int
      lmp_omp_num_per_job is the openmp number for each lammps job.
    lmp_job_per_node: int
      lmp_job_per_node is the number of lammps job in one node.
    proc_num_per_node: 1-d int list
      proc_num_per_node is the number of processors in each node.
    ssh: bool
      ssh is whether to ssh to computational node.
    host: 1-d string list
      host is the name of computational nodes.
  Returns:
    none
  '''

  host_name_proc = []
  for l in range(len(host)):
    host_name_proc.append(''.join((str(proc_num_per_node[l]), '/', host[l])))
  host_info = data_op.comb_list_2_str(host_name_proc, ',')

  run_1 = '''
#! /bin/bash

sys_job="%s"
task_job="%s"
mpi_num="%s"
device_num="%s"
device_start="%s"
direc=%s
parallel_exe=%s

sys_job_arr=(${sys_job///})
task_job_arr=(${task_job///})
mpi_num_arr=(${mpi_num///})
device_num_arr=(${device_num///})
device_start_arr=(${device_start///})

num=${#sys_job_arr[*]}

for ((i=0;i<=num-1;i++));
do
sys_task_mpi_num_arr[i]="${sys_job_arr[i]} ${task_job_arr[i]} ${mpi_num_arr[i]} ${device_num_arr[i]} ${device_start_arr[i]}"
done
''' %(sys_index_str, task_index_str, mpi_num_str, device_num_str, device_id_start_str, lmp_dir, parallel_exe)
  if ssh:
    run_2 = '''
for i in "${sys_task_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --controlmaster -S %s --sshdelay 0.2 $direc/produce.sh {} $direc
''' %(lmp_job_per_node, host_info)
  else:
    run_2 = '''
for i in "${sys_task_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --delay 0.2 $direc/produce.sh {} $direc
''' %(lmp_job_per_node)

  produce = '''
#! /bin/bash

x=$1
direc=$2

x_arr=(${x///})

new_direc=$direc/sys_${x_arr[0]}/task_${x_arr[1]}

lmp_path=%s
mpi_path=%s

ulimit -u 204800
export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

device_num=${x_arr[3]}
if [ $device_num != 0 ]; then
device_id_start=${x_arr[4]}
for ((i=0;i<=device_num-1;i++));
do
((m[i]=$i+device_id_start))
done

for ((i=0;i<device_num;i++));
do
str=$str${m[i]}","
done

export CUDA_VISIBLE_DEVICES=${str:0:(2*$device_num-1)}
fi

a=0

while true
do
if [[ -f $new_direc/atom$a.dump && -f $new_direc/lammps$a.out ]]
then
((a=$a+1))
else
break
fi
done

export OMP_NUM_THREADS=%d

cd $new_direc
mpirun -np ${x_arr[2]} %s < ./md_in.lammps 1> lammps$a.out 2> lammps.err
cd $direc
''' %(lmp_path, mpi_path, lmp_omp_num_per_job, lmp_exe)

  run_file_name_abs = ''.join((lmp_dir, '/run.sh'))
  with open(run_file_name_abs, 'w') as f:
    f.write(run_1+run_2)

  produce_file_name_abs = ''.join((lmp_dir, '/produce.sh'))
  with open(produce_file_name_abs, 'w') as f:
    f.write(produce)

  subprocess.run('chmod +x run.sh', cwd=lmp_dir, shell=True)
  subprocess.run('chmod +x produce.sh', cwd=lmp_dir, shell=True)
  try:
    subprocess.run("bash -c './run.sh'", cwd=lmp_dir, shell=True)
  except subprocess.CalledProcessError as err:
    log_info.log_error('Running error: %s command running error in %s' %(err.cmd, lmp_dir))
    exit()

def run_lmpmd(work_dir, iter_id, lmp_path, lmp_exe, parallel_exe, mpi_path, lmp_job_per_node, \
              lmp_mpi_num_per_job, lmp_omp_num_per_job, proc_num_per_node, host, ssh, device):

  '''
  rum_lmpmd: kernel function to run lammps md.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    lmp_path: string
      lmp_path is the path of lammps.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    parallel_exe: string
      parallel_exe is the parallel executable file.
    mpi_path: string
      mpi_path is the path of mpi.
    lmp_job_per_node: int
      lmp_job_per_node is the number of lammps job in one node.
    lmp_mpi_num_per_job: int
      lmp_mpi_num_per_job is the mpi number for each lammps job.
    lmp_omp_num_per_job: int
      lmp_omp_num_per_job is the openmp number for each lammps job.
    proc_num_per_node: 1-d int list
      proc_num_per_node is the number of processors in each node.
    host: 1-d string list
      host is the name of computational nodes.
    ssh: bool
      ssh is whether we need to ssh.
    device: 2-d int list
      device is the name of gpu devices for all nodes.
  Returns:
    none
  '''

  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(lmp_dir, cmd))

  #All systems have same number of tasks
  lmp_sys_0_dir = ''.join((lmp_dir, '/sys_0'))
  cmd = "ls | grep %s" % ('task_')
  task_num = len(call.call_returns_shell(lmp_sys_0_dir, cmd))

  #check generating lammps tasks
  check_lmp_md_gen = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      lmp_in_file_name_abs = ''.join((lmp_sys_task_dir, '/md_in.lammps'))
      if ( os.path.exists(lmp_in_file_name_abs) and os.path.getsize(lmp_in_file_name_abs) != 0 ):
        check_lmp_md_gen.append(0)
      else:
        check_lmp_md_gen.append(1)

  if ( len(check_lmp_md_gen) !=0 and all(i == 0 for i in check_lmp_md_gen) ):
    str_print = 'Success: generate lammps molecular dynamics tasks in %s' %(lmp_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating lammps molecular dynamics tasks error, please check iteration %d' %(iter_id))
    exit()

  #run lammps md

  if ( sys_num == 1 and task_num == 1 ):
    lmpmd_single(lmp_dir, 0, 0, lmp_exe, lmp_path, mpi_path, lmp_mpi_num_per_job, lmp_omp_num_per_job, device)
  else:
    total_task_num = sys_num*task_num
    sys_task_index = []
    for i in range(sys_num):
      for j in range(task_num):
        sys_task_index.append([i,j])
    calculated_id = 0
    for i in range(total_task_num):
      lmp_sys_task_dir = ''.join((lmp_dir, '/sys_', str(sys_task_index[i][0]), '/task_', str(sys_task_index[i][1])))
      log_file_name = ''.join((lmp_sys_task_dir, '/lammps.out'))
      dump_file_name = ''.join((lmp_sys_task_dir, '/atom.dump'))
      if ( os.path.exists(log_file_name) and os.path.exists(dump_file_name) ):
        calculated_id = i
      else:
        break
    device_num_per_node = [len(device_per_node) for device_per_node in device]
    device_num_min = min(device_num_per_node)
    if ( device_num_min == 1 ):
      lmp_job_per_node = 1

    if ( calculated_id != total_task_num-1 ):
      run_start = calculated_id
      run_end = run_start+lmp_job_per_node*len(host)-1
      if ( run_end > total_task_num-1 ):
        run_end=total_task_num-1
      cycle = math.ceil((total_task_num-run_start)/(lmp_job_per_node*len(host)))
      for i in range(cycle):
        device_num_list = data_op.int_split(device_num_min, lmp_job_per_node)
        device_id_start = [0]
        for j in range(len(device_num_list)-1):
          device_id_start.append(device_id_start[j]+device_num_list[j])
        device_num_str = data_op.comb_list_2_str((device_num_list*len(host))[0:(run_end-run_start+1)], ' ')
        device_id_start_str = data_op.comb_list_2_str((device_id_start*len(host))[0:(run_end-run_start+1)], ' ')

        tot_mpi_num_list = []
        for proc_num in proc_num_per_node:
          mpi_num_list = [lmp_mpi_num_per_job]*lmp_job_per_node
          tot_mpi_num_list.append(mpi_num_list)
        tot_mpi_num_list = data_op.list_reshape(tot_mpi_num_list)[0:(run_end-run_start+1)]
        mpi_num_str = data_op.comb_list_2_str(tot_mpi_num_list, ' ')
        sys_task_index_part = sys_task_index[run_start:run_end+1]
        sys_index = [sys_task[0] for sys_task in sys_task_index_part]
        task_index = [sys_task[1] for sys_task in sys_task_index_part]
        sys_index_str = data_op.comb_list_2_str(sys_index, ' ')
        task_index_str = data_op.comb_list_2_str(task_index, ' ')

        lmpmd_parallel(lmp_dir, lmp_path, mpi_path, lmp_exe, parallel_exe, sys_index_str, task_index_str, mpi_num_str, \
                       device_num_str, device_id_start_str, lmp_omp_num_per_job, lmp_job_per_node, proc_num_per_node, ssh, host)
        for j in range(run_start, run_end+1, 1):
          lmp_sys_task_dir_j = ''.join((lmp_dir, '/sys_', str(sys_task_index[j][0]), '/task_', str(sys_task_index[j][1])))
          gen_lammps_task.combine_frag_traj_file(lmp_sys_task_dir_j)

        run_start = run_start + lmp_job_per_node*len(host)
        run_end = run_end + lmp_job_per_node*len(host)
        if ( run_end > total_task_num-1):
          run_end = total_task_num-1
    else:
      lmpmd_single(lmp_dir, sys_task_index[calculated_id][0], sys_task_index[calculated_id][1], \
                   lmp_exe, lmp_path, mpi_path, lmp_mpi_num_per_job, lmp_omp_num_per_job, device)
  #check lammps md
  check_lmp_md_run = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      dump_file_name_abs = ''.join((lmp_sys_task_dir, '/atom.dump'))
      log_file_name_abs = ''.join((lmp_sys_task_dir, '/lammps.out'))
      if ( os.path.exists(dump_file_name_abs) and \
           os.path.getsize(dump_file_name_abs) != 0 and \
           os.path.exists(log_file_name_abs) and \
           file_tools.grep_line_num("'Step'", log_file_name_abs, lmp_sys_task_dir) != 0 and \
           file_tools.grep_line_num("'Loop time'", log_file_name_abs, lmp_sys_task_dir) != 0 ):
        check_lmp_md_run.append(0)
      else:
        check_lmp_md_run.append(1)

  if ( len(check_lmp_md_run) != 0 and  all(i == 0 for i in check_lmp_md_run) ):
    print ('  Success: molecular dynamics calculations for %d systems by lammps' %(sys_num))
  else:
    log_info.log_error('Running error: lammps molecular dynamics error, please check iteration %d' %(iter_id))
    exit()

def lmpfrc_parallel(model_dir, work_dir, task_index, parallel_exe, lmp_path, \
                    lmp_exe, mpi_path, lmp_job_per_node, host_name, ssh):

  '''
  lmpfrc_parallel: run lammps force calculation in parallel.

  Args:
    model_dir: string
      model_dir is directory for any model.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    task_index: 1-d int list
      task_index is the index of tasks.
    parallel_exe: string
      parallel_exe is the parallel exacutable file.
    lmp_path: string
      lmp_path is the path of lammps.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    mpi_path: string
      mpi_path is the path of mpi.
    lmp_job_per_node: int
      lmp_job_per_node is the number lammps job in each node.
    host_name: string
      host_name is the string host name.
    ssh: bool
      ssh is whether we need to ssh.
  Returns:
    none
  '''

  #run lammps in 1 thread. Here we just run force, it is a single
  #point calculation.

  task_index_str = data_op.comb_list_2_str(task_index, ' ')
  run_1 = '''
#! /bin/bash

direc=%s
task_index="%s"
parallel_exe=%s

task_index_arr=(${task_index///})
num=${#task_index_arr[*]}
''' %(model_dir, task_index_str, parallel_exe)

  if ssh:
    run_2 ='''
for i in "${task_index_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --controlmaster -S %s --sshdelay 0.2  $direc/produce.sh {} $direc
''' %(lmp_job_per_node, host_name)
  else:
    run_2 ='''
for i in "${task_index_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --delay 0.2 $direc/produce.sh {} $direc
''' %(lmp_job_per_node)

  produce = '''
#! /bin/bash

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

x=$1
direc=$2
new_direc=$direc/traj_$x
cd $new_direc
while true
do
%s < $new_direc/frc_in.lammps 1> $new_direc/lammps.out 2> $new_direc/lammps.err
grep "Loop time" $new_direc/lammps.out 1> /dev/null 2> /dev/null
if [[ $? -eq 0 && -f $new_direc/atom.dump ]]
then
break
fi
done
cd %s
''' %(lmp_path, mpi_path, lmp_exe, work_dir)

  produce_file_name_abs = ''.join((model_dir, '/produce.sh'))
  with open(produce_file_name_abs, 'w') as f:
    f.write(produce)
  run_file_name_abs = ''.join((model_dir, '/run.sh'))
  with open(run_file_name_abs, 'w') as f:
    f.write(run_1+run_2)

  subprocess.run('chmod +x produce.sh', cwd=model_dir, shell=True)
  subprocess.run('chmod +x run.sh', cwd=model_dir, shell=True)
  try:
    subprocess.run("bash -c './run.sh'", cwd=model_dir, shell=True)
  except subprocess.CalledProcessError as err:
    log_info.log_error('Running error: %s command running error in %s' %(err.cmd, model_dir))

def check_lmpfrc(lmp_dir, sys_num, use_mtd_tot, atoms_num_tot):

  '''
  check_lmpfrc: check the statu of lammps force calculation

  Args:
    lmp_dir: string
      lmp_dir is the directory of lammps calculation.
    sys_num: int
      sys_num is the number of systems.
    use_mtd_tot: bool
      use_myd_tot is whethet using metadynamics for whole systems.
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example: {1:192,2:90}
  Returns:
    check_lmp_frc_run: 4-d int list
      check_lmp_frc_run is the statu of lammps force calculation.
  '''

  check_lmp_frc_run = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" % ('task_')
    task_num = process.get_task_num(lmp_sys_dir)
    check_lmp_frc_run_i = []
    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      model_num = process.get_lmp_model_num(lmp_sys_task_dir)
      check_lmp_frc_run_ij = []
      for k in range(model_num):
        model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))
        traj_num = process.get_traj_num(model_dir)
        check_lmp_frc_run_ijk = []
        for l in range(traj_num):
          traj_dir = ''.join((model_dir, '/traj_', str(l)))
          dump_file_name_abs = ''.join((traj_dir, '/atom.dump'))
          log_file_name_abs = ''.join((traj_dir, '/lammps.out'))
          if ( os.path.exists(dump_file_name_abs) and \
               len(open(dump_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+9 ):
            if ( use_mtd_tot[i] or k != 0 ):
              if ( os.path.exists(log_file_name_abs) and \
                   file_tools.grep_line_num("'Step'", log_file_name_abs, traj_dir) != 0 and \
                   file_tools.grep_line_num("'Loop time'", log_file_name_abs, traj_dir) != 0 ):
                check_lmp_frc_run_ijk.append(0)
              else:
                check_lmp_frc_run_ijk.append(1)
            else:
              check_lmp_frc_run_ijk.append(0)
          else:
            check_lmp_frc_run_ijk.append(1)
        check_lmp_frc_run_ij.append(check_lmp_frc_run_ijk)
      check_lmp_frc_run_i.append(check_lmp_frc_run_ij)
    check_lmp_frc_run.append(check_lmp_frc_run_i)

  return check_lmp_frc_run

def run_lmpfrc(work_dir, iter_id, lmp_path, lmp_exe, mpi_path, parallel_exe, \
               proc_num_per_node, host, ssh, atoms_num_tot, use_mtd_tot, active_type):

  '''
  rum_lmpfrc: kernel function to run lammps force calculation.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    lmp_path: string
      lmp_path is the path of lammps.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    mpi_path: string
      mpi_path is the path of mpi.
    parallel_exe: string
      parallel_exe is parallel exacutable file.
    proc_num_per_node: 1-d int list
      proc_num_per_node is the number of processors in each node.
    host: 1-d string list
      host is the name of computational nodes.
    ssh: bool
      ssh is whether we need to ssh.
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example: {1:192,2:90}
    use_mtd_tot: bool
      use_mtd_tot is whethet using metadynamics for whole systems.
    active_type: string
      active_type is the type of active learning.
  Returns :
    none
  '''

  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = 'ls | grep %s' % ('sys_')
  sys_num = len(call.call_returns_shell(lmp_dir, cmd))

  #Check generating lammps force tasks.
  check_lmp_frc_gen = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
    task_num = process.get_task_num(lmp_sys_dir)
    if ( active_type == 'model_devi' or use_mtd_tot[i] ):
      for j in range(task_num):
        lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
        model_num = process.get_lmp_model_num(lmp_sys_task_dir)
        for k in range(model_num):
          model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))
          traj_num = process.get_traj_num(model_dir)
          for l in range(traj_num):
            traj_dir = ''.join((model_dir, '/traj_', str(l)))
            if ( not use_mtd_tot[i] and k == 0 ) :
              check_file_name_abs = ''.join((traj_dir, '/atom.dump'))
            else:
              check_file_name_abs = ''.join((traj_dir, '/frc_in.lammps'))
            if ( os.path.exists(check_file_name_abs) and os.path.getsize(check_file_name_abs) != 0 ):
              check_lmp_frc_gen.append(0)
            else:
              check_lmp_frc_gen.append(1)
    else:
      for j in range(task_num):
        check_lmp_frc_gen.append(0)

  if ( len(check_lmp_frc_gen) != 0 and all(i == 0 for i in check_lmp_frc_gen) ):
    str_print = 'Success: generating lammps force calculation file in %s' %(lmp_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating lammps force calculation tasks error, please check iteration %d' %(iter_id))
    exit()

  #Run lammps force.
  for i in range(sys_num):
    if ( active_type == 'model_devi' or use_mtd_tot[i] ):
      lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
      task_num = process.get_task_num(lmp_sys_dir)
      for j in range(task_num):
        lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
        model_num = process.get_lmp_model_num(lmp_sys_task_dir)
        for k in range(model_num):
          model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))
          host_name_proc = []
          for l in range(len(host)):
            host_name_proc.append(''.join((str(proc_num_per_node[l]), '/', host[l])))
          host_info = data_op.comb_list_2_str(host_name_proc, ',')
          traj_num = process.get_traj_num(model_dir)

          calculated_id = 0
          undo_task = []
          for l in range(traj_num):
            traj_dir = ''.join((model_dir, '/traj_', str(l)))
            dump_file_name_abs = ''.join((traj_dir, '/atom.dump'))
            log_file_name_abs = ''.join((traj_dir, '/lammps.out'))
            if ( os.path.exists(dump_file_name_abs) and \
                 len(open(dump_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+9 and \
                 os.path.exists(log_file_name_abs) and \
                 file_tools.grep_line_num("'Step'", log_file_name_abs, traj_dir) != 0 and \
                 file_tools.grep_line_num("'Loop time'", log_file_name_abs, traj_dir) != 0 ):
              pass
            else:
              undo_task.append(l)

          start = 0
          end = min(len(undo_task), start+int(sum(proc_num_per_node)/2))
          cycle = math.ceil(len(undo_task)/int(sum(proc_num_per_node)/2))

          for l in range(cycle):
            if ( use_mtd_tot[i] or k != 0 ):
              task_index = undo_task[start:end]
              lmpfrc_parallel(model_dir, work_dir, task_index, parallel_exe, lmp_path, lmp_exe, \
                                mpi_path, int(proc_num_per_node[0]/2), host_info, ssh)
            start = min(len(undo_task)-1, start+int(sum(proc_num_per_node)/2))
            end = min(len(undo_task), end+int(sum(proc_num_per_node)/2))

  for cycle_run in range(100):
    check_lmp_frc_run = check_lmpfrc(lmp_dir, sys_num, use_mtd_tot, atoms_num_tot)
    lmp_frc_statu = []
    for i in range(sys_num):
      lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
      task_num = process.get_task_num(lmp_sys_dir)
      if ( active_type == 'model_devi' or use_mtd_tot[i] ):
        for j in range(task_num):
          for k in range(model_num):
            if ( len(check_lmp_frc_run[i][j][k]) != 0 and all(m == 0 for m in check_lmp_frc_run[i][j][k]) ):
              lmp_frc_statu.append(0)
            else:
              lmp_frc_statu.append(1)
              failure_task_id = [index for (index,value) in enumerate(check_lmp_frc_run[i][j][k]) if value==1]
              model_dir = ''.join((lmp_dir, '/sys_', str(i), '/task_', str(j), '/model_', str(k)))
              cycle = math.ceil(len(failure_task_id)/int(sum(proc_num_per_node)/2))
              start = 0
              end = min(len(failure_task_id), int(sum(proc_num_per_node)/2))
              for l in range(cycle):
                task_index = failure_task_id[start:end]
                lmpfrc_parallel(model_dir, work_dir, task_index, parallel_exe, lmp_path, lmp_exe, \
                                mpi_path, int(proc_num_per_node[0]/2), host_info, ssh)
                start = min(len(failure_task_id)-1, start+int(sum(proc_num_per_node)/2))
                end = min(len(failure_task_id), end+int(sum(proc_num_per_node)/2))
      else:
        for j in range(task_num):
          lmp_frc_statu.append(0)
    if ( len(lmp_frc_statu) !=0 and all(i == 0 for i in lmp_frc_statu) ):
      break
    if ( cycle_run == 99 and not all(i == 0 for i in lmp_frc_statu) ):
      lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
      task_num = process.get_task_num(lmp_sys_dir)
      for j in range(task_num):
        for k in range(model_num):
          failure_task_id = [index for (index,value) in enumerate(check_lmp_frc_run[i][j][k]) if value==1]
          if ( len(failure_task_id) != 0 ):
            failure_task_id_str = data_op.comb_list_2_str(failure_task_id, ' ')
            str_print = '  Warning: force calculations fail for tasks %s in system %d in task %d by model %d' %(failure_task_id_str, i, j, k)
            str_print = data_op.str_wrap(str_print, 80, '  ')
            print (str_print, flush=True)
      exit()

  print ('  Success: lammps force calculations for %d systems by lammps' %(sys_num), flush=True)

if __name__ == '__main__':
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import lammps_run
  from CP2K_kit.deepff import check_deepff

  box_file = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md/lmp_init_data/box'
  coord_file = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md/lmp_init_data/str.inc'
  tri_cell_vec, atoms, x, y, z = get_box_coord(box_file, coord_file)
  print (atoms, x, y, z)

  exit()
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  deepmd_dic, lmp_dic, cp2k_dic, model_devi_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  proc_num = 4
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  #Test gen_lmpmd_task function
  atoms_type_multi_sys, atoms_num_tot = \
  lmp_run.gen_lmpmd_task(lmp_dic, work_dir, 0)
  #print (atoms_type_multi_sys, atoms_num_tot)

  #Test run_lmpmd function
  lmp_run.run_lmpmd(work_dir, 0, 4)

  #Test gen_lmpfrc_file function
  atoms_num_tot = {0:3}
  atoms_type_multi_sys = {0: {'O': 1, 'H': 2}}
  lmp_run.gen_lmpfrc_file(work_dir, 0, atoms_num_tot, atoms_type_multi_sys)

  #Test run_lmpfrc
  lmp_run.run_lmpfrc(work_dir, 0, 4)
