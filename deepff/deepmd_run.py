#! /usr/env/bin python

import os
import math
import json
import subprocess
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op

def deepmd_parallel_cycle(work_dir, deepmd_train_dir, dp_path, cuda_dir, dp_cmd, dp_version, model_str, \
                          device_num_str, device_start_str, dp_job_per_node, parallel_exe, host):

  '''
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    deepmd_train_dir: string
      deepmd_train_dir is the directory of deepmd-kit training.
    dp_path: string
      dp_path is the path of deepmd-kit.
    cuda_dir: string
      cuda_dir is the directory of cuda.
    dp_cmd: string
      dp_cmd is the deepmd-kit command.
    dp_version: string
      dp_version is the version of deepmd-kit.
    model_str: string
      model_str is the string containing model index.
    device_num_str: string
      device_num_str is the string containing gpu devices number.
    device_id_start_str: string
      device_id_start_str is the string containing staring gpu device id.
    dp_job_per_node: int
      dp_job_per_node is the number of deepmd-kit jobs in one node.
    parallel_exe: string
      parallel_exe is the parallel exacutable file.
    host: 1-d string list
      host is the name of computational nodes.
  '''

  host_info = data_op.comb_list_2_str(host, ',')

  run_1 = '''
#! /bin/bash

model="%s"
device_num="%s"
device_start="%s"
direc=%s
parallel_exe=%s

model_arr=(${model///})
device_num_arr=(${device_num///})
device_start_arr=(${device_start///})

num=${#model_arr[*]}

for ((i=0;i<=num-1;i++));
do
model_device_arr[i]="${model_arr[i]} ${device_num_arr[i]} ${device_start_arr[i]}"
done
''' %(model_str, device_num_str, device_start_str, deepmd_train_dir, parallel_exe)

  if ssh:
    run_2 = '''
for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --controlmaster -S %s --sshdelay 0.2 $direc/produce.sh {} $direc
''' %(dp_job_per_node, host_info)
  else:
    run_2 = '''
for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --delay 0.2 $direc/produce.sh {} $direc
''' %(dp_job_per_node)

  produce_1 = '''
#! /bin/bash

dp_path=%s
export PATH=$dp_path/bin:$PATH
export LD_LIBRARY_PATH=$dp_path/lib:$LD_LIBRARY_PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

x=$1
direc=$2

x_arr=(${x///})

device_num=${x_arr[1]}
if [ $device_num != 0 ]; then
device_id_start=${x_arr[2]}
for ((i=0;i<=device_num-1;i++));
do
((m[i]=$i+device_id_start))
done

for ((i=0;i<device_num;i++));
do
str=$str${m[i]}","
done
''' %(dp_path, cuda_dir)

  if ( dp_version == '1.3.3' ):
    produce_2 ='''
cd $direc/${x_arr[0]}
CUDA_VISIBLE_DEVICES=${str:0:(2*$device_num-1)} \
%s
dp freeze -o frozen_model.pb 1>> log.err 2>> log.err
cd %s
''' %(dp_cmd, work_dir)
  elif ( dp_version == '2.0.0' ):
    produce_2 ='''
cd $direc/${x_arr[0]}
CUDA_VISIBLE_DEVICES=${str:0:(2*$device_num-1)} horovodrun -np $device_num \
%s
dp freeze -o frozen_model.pb 1>> log.err 2>> log.err
cd %s
''' %(dp_cmd, work_dir)

  run_file_name_abs = ''.join((deepmd_train_dir, '/run.sh'))
  with open(run_file_name_abs, 'w') as f:
    f.write(run_1+run_2)

  produce_file_name_abs = ''.join((deepmd_train_dir, '/produce.sh'))
  with open(produce_file_name_abs, 'w') as f:
    f.write(produce_1+produce_2)

  subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
  subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
  try:
    subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
  except subprocess.CalledProcessError as err:
    log_info.log_error('Running error: %s command running error in %s' %(err.cmd, deepmd_train_dir))
    exit()

def deepmd_parallel(work_dir, iter_id, use_prev_model, start, end, parallel_exe, dp_path, host, device, cuda_dir, dp_version, dp_job_per_node):

  '''
  deepmd_parallel: run deepmd calculation in parallel.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    use_prev_model: bool
      use_prev_model is whether we need to use previous model.
    start: int
      start is the starting model id.
    end: int
      end is the endding model id.
    parallel_exe: string
      parallel_exe is the parallel exacutable file.
    dp_path: string
      dp_path is the path of deepmd-kit.
    host: 1-d string list
      host is the name of computational nodes.
    device: 2-d string list
      device is the GPU device.
    cuda_dir: string
      cuda_dir is the directory of cuda.
    dp_version: string
      dp_version is the version of deepmd-kit.
    dp_job_per_node: int
      dp_job_per_node is the number of deepmd-kit jobs in one node.
  Returns:
    none
  '''

  #run lammps in 1 thread. Here we just run force, it is a single
  #point calculation.

  model_num = end-start+1

  deepmd_train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  model_ckpt_file_exists = []
  lcurve_file_exists = []
  final_batch = []
  for i in range(model_num):
    if ( i == 0 ):
      input_file = ''.join((deepmd_train_dir, '/', str(i), '/input.json'))
      with open(input_file, 'r') as f:
        deepmd_dic = json.load(f)
      save_freq = deepmd_dic['training']['save_freq']
    model_ckpt_file = ''.join((deepmd_train_dir, '/', str(i), '/model.ckpt.index'))
    lcurve_file = ''.join((deepmd_train_dir, '/', str(i), '/lcurve.out'))
    model_ckpt_file_exists.append(os.path.exists(model_ckpt_file))
    lcurve_file_exists.append(os.path.exists(lcurve_file))
    if ( os.path.exists(lcurve_file) ):
      whole_line_num = len(open(lcurve_file, 'r').readlines())
      line = linecache.getline(lcurve_file, whole_line_num)
      linecache.clearcache()
      line_split = data_op.split_str(line, ' ', '\n')
      if ( len(line_split) > 0 and data_op.eval_str(line_split[0]) == 1 ):
        final_batch.append(int(line_split[0]))
      else:
        final_batch.append(0)
    else:
      final_batch.append(0)

  if ( all(file_exist for file_exist in model_ckpt_file_exists) and \
       all(file_exist for file_exist in lcurve_file_exists) and \
       all(i > save_freq for i in final_batch) ):
    if ( dp_version == '1.3.3' ):
      dp_cmd = 'dp train --restart model.ckpt input.json 1>> log.out 2>> log.err'
    elif ( dp_version == '2.0.0' ):
      dp_cmd = 'dp train --mpi-log=workers --restart model.ckpt input.json 1>> log.out 2>> log.err'
  else:
    if ( use_prev_model and iter_id>0 ):
      for i in range(model_num):
        model_dir = ''.join((deepmd_train_dir, '/', str(i)))
        prev_model_dir = ''.join((work_dir, '/iter_', str(iter_id-1), '/01.train/', str(i)))
        cmd = "cp %s %s" %('model.ckpt.*', model_dir)
        call.call_simple_shell(prev_model_dir, cmd)
      if ( dp_version == '1.3.3' ):
        dp_cmd = 'dp train --init-model model.ckpt input.json 1> log.out 2> log.err'
      elif ( dp_version == '2.0.0' ):
        dp_cmd = 'dp train --mpi-log=workes --init-model model.ckpt input.json 1> log.out 2> log.err'
    else:
      if ( dp_version == '1.3.3' ):
        dp_cmd = 'dp train input.json 1> log.out 2> log.err'
      elif ( dp_version == '2.0.0' ):
        dp_cmd = 'dp train --mpi-log=workers input.json 1> log.out 2> log.err'

  #Case 1
  if ( len(host) == 1 and len(device[0]) == 0 ):
    run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

seq $run_start $run_end | $parallel_exe -j $parallel_num $direc/produce.sh {} $direc
''' %(deepmd_train_dir, model_num, start, end, parallel_exe)

    produce = '''
#! /bin/bash

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

dp_path=%s

export PATH=$dp_path/bin:$PATH
export LD_LIBRARY_PATH=$dp_path/lib:$LD_LIBRARY_PATH

x=$1
direc=$2
cd $direc/$x
%s
dp freeze -o frozen_model.pb 1>> log.err 2>> log.err
cd %s
''' %(dp_path, dp_cmd, work_dir)

    run_file_name_abs = ''.join((deepmd_train_dir, '/run.sh'))
    with open(run_file_name_abs, 'w') as f:
      f.write(run)

    produce_file_name_abs = ''.join((deepmd_train_dir, '/produce.sh'))
    with open(produce_file_name_abs, 'w') as f:
      f.write(produce)

    subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
    try:
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
    except subprocess.CalledProcessError as err:
      log_info.log_error('Running error: %s command running error in %s' %(err.cmd, deepmd_train_dir))

  #Case 2
  if ( len(host) > 1 and all(len(i) == 0 for i in device) ):
    host_list = []
    for i in range(len(host)):
      host_list.append('-S' + ' ' + host[i])
    if ( len(host) < model_num ):
      host_list = host_list*math.ceil(model_num/len(host))
    host_list = host_list[0:model_num]
    host_comb = data_op.comb_list_2_str(host_list, ' ')

    run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

seq $run_start $run_end | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(deepmd_train_dir, math.ceil(model_num/len(host)), start, end, parallel_exe, host_comb)

    produce = '''
#! /bin/bash

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

dp_path=%s

export PATH=$dp_path/bin:$PATH
export LD_LIBRARY_PATH=$dp_path/lib:$LD_LIBRARY_PATH

x=$1
direc=$2
cd $direc/$x
%s
dp freeze -o frozen_model.pb 1>> log.err 2>> log.err
cd %s
''' %(dp_path, dp_cmd, work_dir)

    run_file_name_abs = ''.join((deepmd_train_dir, '/run.sh'))
    with open(run_file_name_abs, 'w') as f:
      f.write(run)

    produce_file_name_abs = ''.join((deepmd_train_dir, '/produce.sh'))
    with open(produce_file_name_abs, 'w') as f:
      f.write(produce)

    subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
    try:
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
    except subprocess.CalledProcessError as err:
      log_info.log_error('Running error: %s command running error in %s' %(err.cmd, deepmd_train_dir))

  #Case 3
  #If there is just 1 gpu device in the node, we prefer to run deepmd in serial.
  if ( len(host) == 1 and len(device[0]) == 1 ):
    run = '''
#! /bin/bash

dp_path=%s

export PATH=$dp_path/bin:$PATH
export LD_LIBRARY_PATH=$dp_path/lib:$LD_LIBRARY_PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

%s
dp freeze -o frozen_model.pb 1>> log.err 2>> log.err
''' %(dp_path, cuda_dir, dp_cmd)
    for i in range(model_num):
      deepmd_train_i_dir = ''.join((deepmd_train_dir, '/', str(i)))
      run_file_name_abs = ''.join((deepmd_train_i_dir, '/run.sh'))
      with open(run_file_name_abs, 'w') as f:
        f.write(run)
      subprocess.run('chmod +x run.sh', cwd=deepmd_train_i_dir, shell=True)
      try:
        subprocess.run("bash -c './run.sh'", cwd=deepmd_train_i_dir, shell=True)
      except subprocess.CalledProcessError as err:
        log_info.log_error('Running error: %s command running error in %s' %(err.cmd, deepmd_train_i_dir))

  #Case 4
  if ( len(device[0]) > 1 ):
    device_num_per_node = [len(device_per_node) for device_per_node in device]
    device_num_min = min(device_num_per_node)
    cycle = math.ceil(model_num/(dp_job_per_node*len(host)))

    run_start = 0
    run_end = run_start+dp_job_per_node*len(host)-1
    if ( run_end > model_num-1 ):
      run_end=model_num-1

    for i in range(cycle):
      device_num_list = data_op.int_split(device_num_min, lmp_job_per_node)
      device_id_start = [0]
      for j in range(len(device_num_list)-1):
        device_id_start.append(device_id_start[j]+device_num_list[j])
      device_num_str = data_op.comb_list_2_str((device_num_list*len(host))[0:(run_end-run_start+1)], ' ')
      device_id_start_str = data_op.comb_list_2_str((device_id_start*len(host))[0:(run_end-run_start+1)], ' ')
      model_list = data_op.gen_list(run_start, run_end, 1)
      model_str = data_op.comb_list_2_str(model_list, ' ')
      deepmd_parallel_cycle(work_dir, deepmd_train_dir, dp_path, cuda_dir, dp_cmd, dp_version, model_str, \
                            device_num_str, device_start_id_str, dp_job_per_node, parallel_exe, host)

      run_start = run_start + dp_job_per_node*len(host)
      run_end = run_end + dp_job_per_node*len(host)
      if ( run_end > model_num-1):
        run_end = model_num-1

      run_start = 0
      run_end = run_start+len(device[0])-1

def run_deepmd(work_dir, iter_id, use_prev_model, parallel_exe, dp_path, host, device, cuda_dir, dp_version, dp_job_per_node):

  '''
  run_deepmd: kernel function to run deepmd.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    use_prev_model: bool
      use_prev_model is whether we need to use previous model.
    parallel_exe: string
      parallel_exe is the parallel exacutable file.
    dp_path: string
      dp_path is the path of deepmd-kit.
    host: 1-d string list
      host is the name of computational nodes.
    device: 2-d string list
      device is the GPU device.
    cuda_dir: string
      cuda_dir is the directory of cuda.
    dp_version: string
      dp_version is the version of deepmd-kit.
    dp_job_per_node: int
      dp_job_per_node is the number of deepmd-kit jobs in one node.
  Returns:
    none
  '''

  import subprocess

  train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  model_num = len(call.call_returns_shell(train_dir, "ls -ll |awk '/^d/ {print $NF}'"))

  #Check generating deepmd tasks
  check_deepmd_gen = []
  for i in range(model_num):
    inp_file_name_abs = ''.join((train_dir, '/', str(i), '/input.json'))
    if ( os.path.exists(inp_file_name_abs) and os.path.getsize(inp_file_name_abs) != 0 ):
      check_deepmd_gen.append(0)
    else:
      check_deepmd_gen.append(1)

  if ( len(check_deepmd_gen) != 0 and all(i == 0 for i in check_deepmd_gen) ):
    str_print = 'Success: generate deepmd-kit tasks in %s' %(train_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating deepmd-kit tasks error, please check iteration %d' %(iter_id))
    exit()

  if ( not all(len(i) == 0 for i in device) and cuda_dir == 'none' ):
    log_info.log_error('Input error: there are gpu devices in nodes, but cuda_dir is none, please set cuda directory in deepff/environ/cuda_dir')
    exit()
  #Run deepmd-kit tasks
  deepmd_parallel(work_dir, iter_id, use_prev_model, 0, model_num-1, parallel_exe, dp_path, host, device, cuda_dir, dp_version, dp_job_per_node)

  #Check the deepmd tasks.
  check_deepmd_run = []
  for i in range(model_num):
    ff_file_name_abs = ''.join((train_dir, '/', str(i), '/frozen_model.pb'))
    if ( os.path.exists(ff_file_name_abs) and os.path.getsize(ff_file_name_abs) ):
      check_deepmd_run.append(0)
    else:
      check_deepmd_run.append(1)

  if ( len(check_deepmd_run) != 0 and all(i == 0 for i in check_deepmd_run) ):
    print ('  Success: train %d models by deepmd-kit' %(model_num), flush=True)
  else:
    log_info.log_error('deepmd-kit running error, please check iteration %d' %(iter_id))
    exit()

if __name__ == '__main__':
  from CP2K_kit.deepff import deepmd_run
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import load_data
  from CP2K_kit.deepff import check_deepff

  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'

  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  proc_num = 4
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  seed = [1,2,3,4]
  numb_test = int(deepmd_dic['training']['numb_test'])

  i = 0
  init_train_data = []
  cmd = "mkdir %s" % ('init_train_data')
  call.call_simple_shell(work_dir, cmd)
  train_dic = deepmd_dic['training']
  for key in train_dic:
    if ( 'system' in key):
      init_train_key_dir = train_dic[key]['directory']
      proj_name = train_dic[key]['proj_name']
      save_dir = ''.join((work_dir, '/init_train_data/data_', str(i)))
      init_train_data.append(save_dir)
      load_data.load_data_from_dir(init_train_key_dir, work_dir, save_dir, proj_name)
      load_data.raw_to_set(save_dir, 1)
      i = i+1

  #Test gen_deepmd_task function
  deepmd_run.gen_deepmd_task(deepmd_dic, work_dir, 0, init_train_data, seed, numb_test)
  #Test run_deepmd function
  deepmd_run.run_deepmd(work_dir, 0)
