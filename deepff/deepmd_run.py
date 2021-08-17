#! /usr/env/bin python

import math
import json
import copy
import subprocess
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import list_dic_op

def revise_deepmd_dic(deepmd_param):

  '''
  revise_deepmd_dic : transfrom the deepmd keyword in a proper format.

  Args :
    deepmd_param : dictionary
      deepmd_param contains keywords used in deepmd.
  Returns :
    new_deepmd_param : dictionary
      new_deepmd_param is the new deepmd_param.
  '''

  new_deepmd_param = copy.deepcopy(deepmd_param)

  ##revise model descriptor part
  #####################################################################################
  if ( 'sel' in new_deepmd_param['model']['descriptor'].keys() ):
    sel = new_deepmd_param['model']['descriptor']['sel']
    new_deepmd_param['model']['descriptor']['sel'] = [int(x) for x in sel]
  else:
    print ('Cannot find sel, please set sel in model')
    exit()

  if ( 'rcut_smth' in new_deepmd_param['model']['descriptor'].keys() ):
    rcut_smth = new_deepmd_param['model']['descriptor']['rcut_smth']
    new_deepmd_param['model']['descriptor']['rcut_smth'] = float(rcut_smth)
  else:
    new_deepmd_param['model']['descriptor']['rcut_smth'] = 5.0

  if ( 'rcut' in new_deepmd_param['model']['descriptor'].keys() ):
    rcut = new_deepmd_param['model']['descriptor']['rcut']
    new_deepmd_param['model']['descriptor']['rcut'] = float(rcut)
  else:
    new_deepmd_param['model']['descriptor']['rcut'] = 6.0

  if ( 'neuron' in new_deepmd_param['model']['descriptor'].keys() ):
    neuron_encode = new_deepmd_param['model']['descriptor']['neuron']
    new_deepmd_param['model']['descriptor']['neuron'] = [int(x) for x in neuron_encode]
  else:
    new_deepmd_param['model']['descriptor']['neuron'] = [25, 50, 100]

  if ( 'resnet_dt' in new_deepmd_param['model']['descriptor'].keys() ):
    resnet_dt = new_deepmd_param['model']['descriptor']['resnet_dt']
    resnet_dt_bool = list_dic_op.str_to_bool(resnet_dt)
    new_deepmd_param['model']['descriptor']['resnet_dt'] = resnet_dt_bool
  else:
    new_deepmd_param['model']['descriptor']['resnet_dt'] = False

  if ( 'axis_neuron' in new_deepmd_param['model']['descriptor'].keys() ):
    axis_neuron =new_deepmd_param['model']['descriptor']['axis_neuron']
    new_deepmd_param['model']['descriptor']['axis_neuron'] = int(axis_neuron)
  else:
    new_deepmd_param['model']['descriptor']['axis_neuron'] = 16

  ##################################################################################### 

  ##revise model fitting_net part
  #####################################################################################
  if ( 'resnet_dt' in new_deepmd_param['model']['fitting_net'].keys() ):
    resnet_dt = new_deepmd_param['model']['fitting_net']['resnet_dt']
    resnet_dt_bool = list_dic_op.str_to_bool(resnet_dt)
    new_deepmd_param['model']['fitting_net']['resnet_dt'] = resnet_dt_bool
  else:
    new_deepmd_param['model']['fitting_net']['resnet_dt'] = True

  ##revise learning_rate part
  #####################################################################################
  if ( 'decay_steps' in new_deepmd_param['learning_rate'].keys() ):
    decay_steps = new_deepmd_param['learning_rate']['decay_steps']
    new_deepmd_param['learning_rate']['decay_steps'] = int(decay_steps)
  else:
    new_deepmd_param['learning_rate']['decay_steps'] = 5000

  if ( 'start_lr' in new_deepmd_param['learning_rate'].keys() ):
    start_lr = new_deepmd_param['learning_rate']['start_lr']
    new_deepmd_param['learning_rate']['start_lr'] = float(start_lr)
  else:
    new_deepmd_param['learning_rate']['start_lr'] = 0.001

  if ( 'stop_lr' in new_deepmd_param['learning_rate'].keys() ):
    stop_lr = new_deepmd_param['learning_rate']['stop_lr']
    new_deepmd_param['learning_rate']['stop_lr'] = float(stop_lr)
  else:
    new_deepmd_param['learning_rate']['stop_lr'] = 1e-8
  #####################################################################################

  ##revise loss part
  #####################################################################################
  if ( 'start_pref_e' in new_deepmd_param['loss'].keys() ):
    start_pref_e = new_deepmd_param['loss']['start_pref_e']
    new_deepmd_param['loss']['start_pref_e'] = float(start_pref_e)
  else:
    new_deepmd_param['loss']['start_pref_e'] = 0.02

  if ( 'limit_pref_e' in new_deepmd_param['loss'].keys() ):
    limit_pref_e = new_deepmd_param['loss']['limit_pref_e']
    new_deepmd_param['loss']['limit_pref_e'] = float(limit_pref_e)
  else:
    new_deepmd_param['loss']['limit_pref_e'] = 1.0

  if ( 'start_pref_f' in new_deepmd_param['loss'].keys() ):
    start_pref_f = new_deepmd_param['loss']['start_pref_f']
    new_deepmd_param['loss']['start_pref_f'] = float(start_pref_f)
  else:
    new_deepmd_param['loss']['start_pref_f'] = 1000.0

  if ( 'limit_pref_f' in new_deepmd_param['loss'].keys() ):
    limit_pref_f = new_deepmd_param['loss']['limit_pref_f']
    new_deepmd_param['loss']['limit_pref_f'] = float(limit_pref_f)
  else:
    new_deepmd_param['loss']['limit_pref_f'] = 1.0

  if ( 'start_pref_v' in new_deepmd_param['loss'].keys() ):
    start_pref_v = new_deepmd_param['loss']['start_pref_v']
    new_deepmd_param['loss']['start_pref_v'] = float(start_pref_v)
  else:
    new_deepmd_param['loss']['start_pref_v'] = 0.0

  if ( 'limit_pref_v' in new_deepmd_param['loss'].keys() ):
    limit_pref_v = new_deepmd_param['loss']['limit_pref_v']
    new_deepmd_param['loss']['limit_pref_v'] = float(limit_pref_v)
  else:
    new_deepmd_param['loss']['limit_pref_v'] = 0.0
  ######################################################################################

  ##revise training part
  ######################################################################################
  if ( 'stop_batch' in new_deepmd_param['training'].keys() ):
    stop_batch = new_deepmd_param['training']['stop_batch']
    new_deepmd_param['training']['stop_batch'] = int(stop_batch)
  else:
    new_deepmd_param['training']['stop_batch'] = 40000

  if ( 'batch_size' in new_deepmd_param['training'].keys() ):
    batch_size = new_deepmd_param['training']['batch_size']
    if ( batch_size.isdigit() ):
      new_deepmd_param['training']['batch_size'] = int(batch_size)
  else:
    new_deepmd_param['training']['batch_size'] = 1

  if ( 'disp_freq' in new_deepmd_param['training'].keys() ):
    disp_freq = new_deepmd_param['training']['disp_freq']
    new_deepmd_param['training']['disp_freq'] = int(disp_freq)
  else:
    new_deepmd_param['training']['disp_freq'] = 100

  if ( 'numb_test' in new_deepmd_param['training'].keys() ):
    numb_test = new_deepmd_param['training']['numb_test']
    new_deepmd_param['training']['numb_test'] = int(numb_test)
  else:
    new_deepmd_param['training']['numb_test'] = 5

  if ( 'save_freq' in new_deepmd_param['training'].keys() ):
    save_freq = new_deepmd_param['training']['save_freq']
    new_deepmd_param['training']['save_freq'] = int(save_freq)
  else:
    new_deepmd_param['training']['save_freq'] = 1000

  if ( 'disp_training' in new_deepmd_param['training'].keys() ):
    disp_training = new_deepmd_param['training']['disp_training']
    disp_training_bool = list_dic_op.str_to_bool(disp_training)
    new_deepmd_param['training']['disp_training'] = disp_training_bool
  else:
    new_deepmd_param['training']['disp_training'] = True

  if ( 'time_training' in new_deepmd_param['training'].keys() ):
    time_training = new_deepmd_param['training']['time_training']
    time_training_bool = list_dic_op.str_to_bool(time_training)
    new_deepmd_param['training']['time_training'] = time_training_bool
  else:
    new_deepmd_param['training']['time_training'] = True

  if ( 'profiling' in new_deepmd_param['training'].keys() ):
    profiling = new_deepmd_param['training']['profiling']
    profiling_bool = list_dic_op.str_to_bool(profiling)
    new_deepmd_param['training']['profiling'] = profiling_bool
  else:
    new_deepmd_param['training']['profiling'] = False
  #####################################################################################

  return new_deepmd_param

def gen_deepmd_task(deepmd_dic, work_dir, iter_id, init_train_data, numb_test, \
                    descr_seed, fit_seed, tra_seed, neuron, model_type):

  '''
  gen_deepmd_task : generate deepmd tasks

  Args :
    deepmd_dic : dictionary
      deepmd_dic contains keywords used in deepmd.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    iter_id : int
      iter_id is the iteration id.
    init_train_data : 1-d string list
      init_train_data contains initial training data directories.
    numb_test : int
      numb_test is number of testing data sets when run deepmd
    seed : 1-d int list or int
      seed is the seed number to generate random number.
    neuron : 1-d int list or int
      neuron is the number of nodes in neural network.
    model_type : string
      model_type has two choices: 'use_seed', 'use_nodes'
  Returns :
    none
  '''

  #copy should be done at first, because the following operators will change it!
  deepmd_param = copy.deepcopy(deepmd_dic)

  iter_dir = ''.join((work_dir, '/iter_', str(iter_id)))
  cmd = "mkdir %s" % ('01.train')
  call.call_simple_shell(iter_dir, cmd)
  train_dir = ''.join(iter_dir + '/01.train')

  data_dir = copy.deepcopy(init_train_data)
  if ( iter_id > 0 ):
    for i in range(iter_id):
      calc_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc'))
      cmd = "ls | grep %s" % ('sys_')
      sys_num = len(call.call_returns_shell(calc_dir, cmd))
      for j in range(sys_num):
        cp2k_sys_dir = ''.join((calc_dir, '/sys_', str(j)))
        cmd = "ls | grep %s" % ('task_')
        task_num = len(call.call_returns_shell(cp2k_sys_dir, cmd))
        if ( task_num > numb_test ):
          prev_data_dir = ''.join((cp2k_sys_dir, '/data'))
          data_dir.append(prev_data_dir)

  for key in deepmd_dic['training']:
    if ( 'system' in key ):
      deepmd_param['training'].pop(key)

  if ( 'model_type' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('model_type')

  if ( 'neuron' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('neuron')

  if ( 'seed_num' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('seed_num')

  deepmd_param['training']['systems'] = data_dir

  new_deepmd_param = revise_deepmd_dic(deepmd_param)

  if ( model_type == 'use_seed' ):
    for i in range(len(descr_seed)):
      cmd = "mkdir %s" % (str(i))
      call.call_simple_shell(train_dir, cmd)
      deepmd_param_i = copy.deepcopy(new_deepmd_param)
      deepmd_param_i['model']['descriptor']['seed'] = descr_seed[i]
      deepmd_param_i['model']['fitting_net']['seed'] = fit_seed[i]
      deepmd_param_i['training']['seed'] = tra_seed[i]
      deepmd_param_i['model']['fitting_net']['neuron'] = neuron
      json_str = json.dumps(deepmd_param_i, indent=4)

      with open(''.join((train_dir, '/', str(i), '/', 'input.json')), 'w') as json_file:
        json_file.write(json_str)

  elif ( model_type == 'use_nodes' ):
    for i in range(len(neuron)):
      cmd = "mkdir %s" % (str(i))
      call.call_simple_shell(train_dir, cmd)
      deepmd_param_i = copy.deepcopy(new_deepmd_param)
      deepmd_param_i['model']['descriptor']['seed'] = descr_seed[i]
      deepmd_param_i['model']['fitting_net']['seed'] = fit_seed[i]
      deepmd_param_i['training']['seed'] = tra_seed[i]
      deepmd_param_i['model']['fitting_net']['neuron'] = neuron[i]
      json_str = json.dumps(deepmd_param_i, indent=4)

      with open(''.join((train_dir, '/', str(i), '/', 'input.json')), 'w') as json_file:
        json_file.write(json_str)

def deepmd_parallel(deepmd_train_dir, start, end, parallel_exe, host, device, usage, cuda_dir):

  '''
  deepmd_parallel : run deepmd calculation in parallel.

  Args :
    deepmd_train_dir : string
      deepmd_train_dir is the directory of deepmd training for each iteration.
    start : int
      start is the starting model id.
    end : int
      end is the endding model id.
    parallel_exe : string
      parallel_exe is the parallel exacutable file.
    host : 1-d string list
      host is the name of computational nodes.
    device : 2-d string list
      device is the GPU device.
    usage : 2-d float list
      usage is the memory use of GPU.
    cuda_dir : cuda_dir is the directory of cuda.
  Returns :
    none
  '''

  #run lammps in 1 thread. Here we just run force, it is a single
  #point calculation.

  model_num = end-start+1
  dp_exe = call.call_returns_shell(deepmd_train_dir, 'which dp')[0]
  dp_path = dp_exe[:-3]

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
export PATH=%s:$PATH

x=$1
direc=$2
cd $direc/$x
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path)

    run_file = ''.join((deepmd_train_dir, '/run.sh'))
    with open(run_file, 'w') as f:
      f.write(run)

    produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
    with open(produce_file, 'w') as f:
      f.write(produce)

    subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)

  #Case 2
  if ( len(host) > 1 and all(len(i) == 0 for i in device) ):
    host_list = []
    for i in range(len(host)):
      host_list.append('-S' + ' ' + host[i])
    if ( len(host) < model_num ):
      host_list = host_list*math.ceil(model_num/len(host))
    host_list = host_list[0:model_num]
    host_comb = list_dic_op.comb_list_2_str(host_list, ' ')

    run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

seq $run_start $run_end | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(deepmd_train_dir, model_num, start, end, parallel_exe, host_comb)

    produce = '''
#! /bin/bash

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PATH=%s:$PATH

x=$1
direc=$2
cd $direc/$x
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path)

    run_file = ''.join((deepmd_train_dir, '/run.sh'))
    with open(run_file, 'w') as f:
      f.write(run)

    produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
    with open(produce_file, 'w') as f:
      f.write(produce)

    subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
    subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)

  #Case 3
  #If there is just 1 gpu device in the node, we prefer to run deepmd in serial.
  if ( len(host) == 1 and len(device[0]) == 1 ):
    run = '''
#! /bin/bash

export PATH=%s:$PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

dp train input.json 1> out 2> err
dp freeze -o frozen_model.pb
''' %(dp_path, cuda_dir)
    for i in range(model_num):
      deepmd_train_i_dir = ''.join((deepmd_train_dir, '/', str(i)))
      run_file = ''.join((deepmd_train_i_dir, '/run.sh'))
      with open(run_file, 'w') as f:
        f.write(run)
      subprocess.run('chmod +x run.sh', cwd=deepmd_train_i_dir, shell=True)
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_i_dir, shell=True)

  #Case 4
  if ( len(host) == 1 and len(device[0]) > 1 ):
    if ( len(device[0]) >= model_num ):
      device_choose = device[0][0:model_num]
      device_str = list_dic_op.comb_list_2_str(device_choose, ' ')
      model_list = list_dic_op.gen_list(start, end, 1)
      model_str = list_dic_op.comb_list_2_str(model_list, ' ')

      run = '''
#! /bin/bash

model="%s"
device="%s"
direc=%s
parallel_num=%d
parallel_exe=%s

model_arr=(${model///})
device_arr=(${device///})

num=${#model_arr[*]}

for ((i=0;i<=num-1;i++));
do
model_device_arr[i]="${model_arr[i]} ${device_arr[i]}"
done

for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j $parallel_num $direc/produce.sh {} $direc
''' %(model_str, device_str, deepmd_train_dir, model_num, parallel_exe)

      produce = '''
#! /bin/bash

export PATH=%s:$PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

x=$1
direc=$2

x_arr=(${x///})

export CUDA_VISIBLE_DEVICES=${x_arr[1]}

cd $direc/${x_arr[0]}
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path, cuda_dir)

      run_file = ''.join((deepmd_train_dir, '/run.sh'))
      with open(run_file, 'w') as f:
        f.write(run)

      produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
      with open(produce_file, 'w') as f:
        f.write(produce)

      subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
    else:
      run_start = 0
      run_end = run_start+len(device[0])-1
      cycle = math.ceil(model_num/len(device[0]))
      for i in range(cycle):
        device_str = list_dic_op.comb_list_2_str(device[0][0:(run_end-run_start+1)], ' ')
        model_list = list_dic_op.gen_list(run_start, run_end, 1)
        model_str = list_dic_op.comb_list_2_str(model_list, ' ')
        run = '''
#! /bin/bash

model="%s"
device="%s"
direc=%s
parallel_num=%d
parallel_exe=%s

model_arr=(${model///})
device_arr=(${device///})

num=${#model_arr[*]}

for ((i=0;i<=num-1;i++));
do
model_device_arr[i]="${model_arr[i]} ${device_arr[i]}"
done

for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j $parallel_num $direc/produce.sh {} $direc
''' %(model_str, device_str, deepmd_train_dir, run_end-run_start+1, parallel_exe)

        produce = '''
#! /bin/bash

export PATH=%s:$PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

x=$1
direc=$2

x_arr=(${x///})

export CUDA_VISIBLE_DEVICES=${x_arr[1]}

cd $direc/${x_arr[0]}
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path, cuda_dir)

        run_file = ''.join((deepmd_train_dir, '/run.sh'))
        with open(run_file, 'w') as f:
          f.write(run)

        produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
        with open(produce_file, 'w') as f:
          f.write(produce)

        subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)

        run_start = run_start + len(device[0])
        run_end = run_end + len(device[0])
        if ( run_end > end):
          run_end = end

  #Case 5
  if ( len(host) > 1 and not all(len(i) == 0 for i in device) and len(list_dic_op.list_replicate(device)) == 1 ):

    host_exp = []
    for i in range(len(host)):
      for j in range(len(device[i])):
        host_exp.append(host[i])
    device_exp = list_dic_op.list_reshape(device)

    if ( len(device_exp) >= model_num ):
      host_list = []
      for i in range(len(model_num)):
        host_list.append('-S' + ' ' + host_exp[i])
      host_list = host_list[0:model_num]
      host_comb = list_dic_op.comb_list_2_str(host_list, ' ')
      device_choose = device_exp[0:model_num]
      device_str = list_dic_op.comb_list_2_str(device_choose, ' ')
      model_list = list_dic_op.gen_list(start, end, 1)
      model_str = list_dic_op.comb_list_2_str(model_list, ' ')

      run = '''
#! /bin/bash

model="%s"
device="%s"
direc=%s
parallel_num=%d
parallel_exe=%s

model_arr=(${model///})
device_arr=(${device///})

num=${#model_arr[*]}

for ((i=0;i<=num-1;i++));
do
model_device_arr[i]="${model_arr[i]} ${device_arr[i]}"
done

for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(model_str, device_str, deepmd_train_dir, model_num, parallel_exe, host_comb)

      produce = '''
#! /bin/bash

export PATH=%s:$PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

x=$1
direc=$2

export CUDA_VISIBLE_DEVICE=${x_arr[1]}

cd $direc/${x_arr[0]}
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path, cuda_dir)

      run_file = ''.join((deepmd_train_dir, '/run.sh'))
      with open(run_file, 'w') as f:
        f.write(run)

      produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
      with open(produce_file, 'w') as f:
        f.write(produce)

      subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
    else:
      run_start = 0
      run_end = run_start+len(device_exp)-1
      cycle = math.ceil(model_num/len(device_exp))
      for i in range(cycle):
        device_str = list_dic_op.comb_list_2_str(device_exp[0:(run_end-run_start+1)], ' ')
        model_list = list_dic_op.gen_list(run_start, run_end, 1)
        model_str = list_dic_op.comb_list_2_str(model_list, ' ')
        host_list = []
        for j in range(run_end-run_start+1):
          host_list.append('-S' + ' ' + host_exp[j])
        host_comb = list_dic_op.comb_list_2_str(host_list, ' ')

        run = '''
#! /bin/bash

model="%s"
device="%s"
direc=%s
parallel_num=%d
parallel_exe=%s

model_arr=(${model///})
device_arr=(${device///})

num=${#model_arr[*]}

for ((i=0;i<=num-1;i++));
do
model_device_arr[i]="${model_arr[i]} ${device_arr[i]}"
done

for i in "${model_device_arr[@]}"; do echo "$i"; done | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(model_str, device_str, deepmd_train_dir, run_end-run_start+1, parallel_exe, host_comb)

        produce = '''
#! /bin/bash

export PATH=%s:$PATH

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=fine,verbose,compact,1,0

x=$1
direc=$2

x_arr=(${x///})

export CUDA_VISIBLE_DEVICE=${x_arr[1]}

cd $direc/${x_arr[0]}
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(dp_path, cuda_dir)

        run_file = ''.join((deepmd_train_dir, '/run.sh'))
        with open(run_file, 'w') as f:
          f.write(run)

        produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
        with open(produce_file, 'w') as f:
          f.write(produce)

        subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)

        run_start = run_start + len(device[0])
        run_end = run_end + len(device[0])
        if ( run_end > end):
          run_end = end

  #Case 6
  if ( len(host) > 1 and not all(len(i) == 0 for i in device) and len(list_dic_op.list_replicate(device)) != 1 ):
    host_true = []
    for i in range(len(host)):
      if ( len(device[i]) != 0 ):
        host_true.append(host[i])
    if ( len(host_true) > model_num ):
      host_list = []
      for i in range(model_num):
        host_list.append('-S' + ' ' + host_true[j])
      host_comb = list_dic_op.comb_list_2_str(host_list, ' ')

      run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

seq $run_start $run_end | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(deepmd_train_dir, model_num, start, end, parallel_exe, host_comb)

      produce = '''
#! /bin/bash

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export PATH=%s:$PATH

x=$1
direc=$2
cd $direc/$x
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(cuda_dir, dp_path)

      run_file = ''.join((deepmd_train_dir, '/run.sh'))
      with open(run_file, 'w') as f:
        f.write(run)

      produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
      with open(produce_file, 'w') as f:
        f.write(produce)

      subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
      subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)
    else:
      run_start = 0
      run_end = run_start+len(host_true)-1
      cycle = math.ceil(model_num/len(host_true))
      for i in range(cycle):
        host_list = []
        for j in range(run_end-run_start+1):
          host_list.append('-S' + ' ' + host_true[j])
        host_comb = list_dic_op.comb_list_2_str(host_list, ' ')
        run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

seq $run_start $run_end | $parallel_exe -j $parallel_num %s $direc/produce.sh {} $direc
''' %(deepmd_train_dir, run_end-run_start+1, run_start, run_end, parallel_exe, host_comb)

        produce = '''
#! /bin/bash

CUDA_DIR=%s
export PATH=$CUDA_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_DIR/lib64:$LD_LIBRARY_PATH

export PATH=%s:$PATH

x=$1
direc=$2
cd $direc/$x
dp train input.json 1> log.out 2> log.err
dp freeze -o frozen_model.pb
''' %(cuda_dir, dp_path)

        run_file = ''.join((deepmd_train_dir, '/run.sh'))
        with open(run_file, 'w') as f:
          f.write(run)

        produce_file = ''.join((deepmd_train_dir, '/produce.sh'))
        with open(produce_file, 'w') as f:
          f.write(produce)

        subprocess.run('chmod +x run.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run('chmod +x produce.sh', cwd=deepmd_train_dir, shell=True)
        subprocess.run("bash -c './run.sh'", cwd=deepmd_train_dir, shell=True)

        run_start = run_start + len(host_true)
        run_end = run_end + len(host_true)
        if ( run_end > end ):
          run_end = end

def run_deepmd(work_dir, iter_id, parallel_exe, host, device, usage, cuda_dir):

  '''
  run_deepmd : kernel function to run deepmd.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    iter_id : int
      iter_id is the iteration id.
    parallel_exe: str
      parallel_exe is the executable file parallel.
  Returns :
    none
  '''

  import subprocess

  train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  model_num = len(call.call_returns_shell(train_dir, "ls -ll |awk '/^d/ {print $NF}'"))

  deepmd_parallel(train_dir, 0, model_num-1, parallel_exe, host, device, usage, cuda_dir)

if __name__ == '__main__':
  from CP2K_kit.deepff import deepmd_run
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import load_data

  deepff_key = ['deepmd', 'lammps', 'cp2k', 'force_eval', 'environ']
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'

  deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
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
