#! /usr/env/bin python

import os
import math
import json
import copy
import numpy as np
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.deepff import process
from CP2K_kit.deepff import write_data

def assign_data_dir(work_dir, init_train_data, iter_id, numb_test):

  '''
  assign_data_dir: assign the directory of data

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    init_train_data: 1-d string list
      init_train_data contains initial training data directories.
    iter_id: int
      iter_id is the iteration id.
    numb_test: int
      numb_test is number of testing data sets when run deepmd.
  Returns:
    data_dir: 1-d string list
      data_dir contains all training data directories.
    final_data_dir: 1-d string list
      final_data_dir contains final iteration training data directories.
  '''

  data_dir = copy.deepcopy(init_train_data)
  final_data_dir = []
  if ( iter_id > 0 ):
    for i in range(iter_id):
      cp2k_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc'))
      sys_num = process.get_sys_num(cp2k_dir)
      for j in range(sys_num):
        cp2k_sys_dir = ''.join((cp2k_dir, '/sys_', str(j)))
        task_num = process.get_task_num(cp2k_sys_dir)
        for k in range(task_num):
          prev_data_dir = ''.join((cp2k_sys_dir, '/task_', str(k), '/data'))
          if ( os.path.exists(prev_data_dir) ):
            cmd = "ls | grep %s" % ('set.')
            set_parts = len(call.call_returns_shell(prev_data_dir, cmd))
            if ( (set_parts-1) >= 10 ):
              final_set_dir_name = 'set.0'+str(set_parts-1)
            elif ( (set_parts-1) < 10 ):
              final_set_dir_name = 'set.00'+str(set_parts-1)
            elif ( (set_parts-1) > 100 ):
              final_set_dir_name = 'set.'+str(set_parts-1)
            ene_npy_file = ''.join((prev_data_dir, '/', final_set_dir_name, '/energy.npy'))
            final_set_data_num = len(np.load(ene_npy_file))
            if ( final_set_data_num > numb_test ):
              data_dir.append(prev_data_dir)
              if ( i == iter_id-1 ):
                final_data_dir.append(prev_data_dir)

  return data_dir, final_data_dir

def assign_prob(data_dir, final_data_dir, data_num, iter_id):

  '''
  assign_prob: assign probability of systems.

  Args:
    data_dir: 1-d string list
      data_dir contains all training data directories.
    final_data_dir: 1-d string list
      final_data_dir contains final iteration training data directories.
    data_num: 1-d int list
      data_num is the numbers of training data for each system.
    iter_id: int
      iter_id is the iteration id.
  Returns:
    prob_sys_1: string
      prob_sys_1 is the probability of system 1.
    prob_sys_2: string
      prob_sys_2 is the probability of system 2.
  '''

  data_num_1 = data_num[0:(len(data_dir)-len(final_data_dir))]
  prob_data_num_1 = math.erf(sum(data_num_1)/sum(data_num)*2)/2**(1/(iter_id+1))
  prob_data_num_2 = 1.0-prob_data_num_1
  prob_sys_1 = '0:%d:%f' %(len(data_dir)-len(final_data_dir), prob_data_num_1)
  prob_sys_2 = '%d:%d:%f' %(len(data_dir)-len(final_data_dir), len(data_dir), prob_data_num_2)
  #prob_sys_1 = '0:%d:%f' %(len(data_dir)-len(final_data_dir), 0.2)
  #prob_sys_2 = '%d:%d:%f' %(len(data_dir)-len(final_data_dir), len(data_dir), 0.8)

  return prob_sys_1, prob_sys_2

def gen_deepmd_model_task(deepmd_dic, work_dir, iter_id, init_train_data, numb_test, descr_seed, fit_seed, \
                          tra_seed, neuron, model_type, data_num, tot_atoms_type_dic):

  '''
  gen_deepmd_model_task: generate deepmd model tasks

  Args :
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    init_train_data: 1-d string list
      init_train_data contains initial training data directories.
    numb_test: int
      numb_test is number of testing data sets when run deepmd
    descr_seed: 1-d int list
      descr_seed is the seed number in descriptor part.
    fit_seed: 1-d int list
      fit_seed is the seed number in fitting part.
    tra_seed: 1-d int list
      tra_seed is the seed number in training part.
    neuron: 1-d int list or int
      neuron is the number of nodes in neural network.
    model_type: string
      model_type has two choices: 'use_seed', 'use_node'
    data_num: 1-d int list
      data_num is the numbers of training data for each system.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  #copy should be done at first, because the following operators will change it!
  deepmd_param = copy.deepcopy(deepmd_dic)

  iter_dir = ''.join((work_dir, '/iter_', str(iter_id)))
  train_dir = ''.join(iter_dir + '/01.train')
  if ( not os.path.exists(train_dir) ):
    cmd = "mkdir %s" % ('01.train')
    call.call_simple_shell(iter_dir, cmd)

  data_dir, final_data_dir = assign_data_dir(work_dir, init_train_data, iter_id, numb_test)

  if ( iter_id > 1 and len(final_data_dir) == 0 ):
    str_print = '''  There are no new data in No.%d iteration, two possible reasons:
  (1) If the accurate ratio is larger than 96.0%%, the deepff is converged
  (2) If the accurate ratio is very small (less than 4.0%%), the deepff is no
  converged. But the program have to be stopped. The whole iterations does not
  succeed. The only thing we can do is changing either initial data sets or
  training parameters.''' %(iter_id-1)
    print ('Warning'.center(80,'*'), flush=True)
    print (str_print, flush=True)
    write_data.write_active_data(work_dir, iter_id, tot_atoms_type_dic)
    exit()

  start_lr = deepmd_param['learning_rate']['start_lr']
  batch_size = deepmd_param['training']['batch_size']
  lr_scale = deepmd_param['training']['lr_scale']
  fix_stop_batch = deepmd_param['training']['fix_stop_batch']
  use_prev_model = deepmd_param['training']['use_prev_model']
  if not fix_stop_batch:
    epoch_num = deepmd_param['training']['epoch_num']
    stop_batch = math.ceil(epoch_num*int(sum(data_num)/batch_size)/10000)*10000
    decay_steps = math.ceil(int(sum(data_num)/batch_size)/1000)*1000
    if ( stop_batch < decay_steps*200 ):
      stop_batch = decay_steps*200
    deepmd_param['learning_rate']['decay_steps'] = decay_steps
    deepmd_param['training']['stop_batch'] = stop_batch

  deepmd_param['training']['systems'] = data_dir
  deepmd_param['training']['batch_size'] = [batch_size]*len(data_dir)

  if ( 'atom_mass' in deepmd_param['model'].keys() ):
    deepmd_param['model'].pop('atom_mass')

  for key in deepmd_dic['training'].keys():
    if ( 'system' in key ):
      deepmd_param['training'].pop(key)

  if ( 'seed_num' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('seed_num')
  if ( 'epoch_num' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('epoch_num')
  if ( 'lr_scale' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('lr_scale')
  if ( 'set_data_dir' in deepmd_param['training'].keys() ):
    deepmd_param['training'].pop('set_data_dir')
  deepmd_param['training'].pop('model_type')
  deepmd_param['training'].pop('neuron')
  deepmd_param['training'].pop('shuffle_data')
  deepmd_param['training'].pop('train_stress')
  deepmd_param['training'].pop('use_prev_model')
  deepmd_param['training'].pop('fix_stop_batch')

  if ( iter_id > 0 and use_prev_model ):
    if ( lr_scale**iter_id > 100.0 ):
      deepmd_param['learning_rate']['start_lr'] = start_lr/(100.0)
    else:
      deepmd_param['learning_rate']['start_lr'] = start_lr/(lr_scale**iter_id)
    deepmd_param['training']['stop_batch'] = int(deepmd_param['training']['stop_batch']/2)
    prob_sys_1, prob_sys_2 = assign_prob(data_dir, final_data_dir, data_num, iter_id)
    auto_prob_style = data_op.comb_list_2_str(['prob_sys_size', prob_sys_1, prob_sys_2], ';')
    deepmd_param['training']['auto_prob_style'] = auto_prob_style

  if ( model_type == 'use_seed' ):
    for i in range(len(descr_seed)):
      cmd = "mkdir %s" % (str(i))
      call.call_simple_shell(train_dir, cmd)
      deepmd_param_i = copy.deepcopy(deepmd_param)
      deepmd_param_i['model']['descriptor']['seed'] = descr_seed[i]
      deepmd_param_i['model']['fitting_net']['seed'] = fit_seed[i]
      deepmd_param_i['training']['seed'] = tra_seed[i]
      deepmd_param_i['model']['fitting_net']['neuron'] = neuron
      json_str = json.dumps(deepmd_param_i, indent=4)

      with open(''.join((train_dir, '/', str(i), '/', 'input.json')), 'w') as json_file:
        json_file.write(json_str)

  elif ( model_type == 'use_node' ):
    for i in range(len(neuron)):
      model_dir = ''.join((train_dir, '/', str(i)))
      if ( not os.path.exists(model_dir) ):
        cmd = "mkdir %s" % (str(i))
        call.call_simple_shell(train_dir, cmd)
      deepmd_param_i = copy.deepcopy(deepmd_param)
      deepmd_param_i['model']['descriptor']['seed'] = descr_seed[i]
      deepmd_param_i['model']['fitting_net']['seed'] = fit_seed[i]
      deepmd_param_i['training']['seed'] = tra_seed[i]
      deepmd_param_i['model']['fitting_net']['neuron'] = neuron[i]
      json_str = json.dumps(deepmd_param_i, indent=4)

      with open(''.join((model_dir, '/input.json')), 'w') as json_file:
        json_file.write(json_str)

def gen_deepmd_test_task(deepmd_test_dic, work_dir, iter_id, data_num, tot_atoms_type_dic):

  '''
  gen_deepmd_test_task: generate deepmd test tasks

  Args :
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    deepmd_inp_file: string
      deepmd_inp_file is the input file of deepmd.
    init_dpff_file: string
      init_dpff_file is the initial deep potential force field file.
    data_num: 1-d int list
      data_num is the numbers of training data for each system.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  iter_dir = ''.join((work_dir, '/iter_', str(iter_id)))
  train_dir = ''.join(iter_dir + '/01.train')
  if ( not os.path.exists(train_dir) ):
    cmd = "mkdir %s" %('01.train')
    call.call_simple_shell(iter_dir, cmd)

  init_dpff_dir = deepmd_test_dic['init_dpff_dir']
  deepmd_inp_file = ''.join((init_dpff_dir, '/input.json'))

  model_dir = ''.join((train_dir, '/0'))
  if ( not os.path.exists(model_dir) ):
    cmd = "mkdir 0"
    call.call_simple_shell(train_dir, cmd)
  if ( iter_id == 0 ):
    cmd = "cp %s %s" %(''.join((init_dpff_dir, '/frozen_model.pb')), model_dir)
    call.call_simple_shell(train_dir, cmd)
    cmd = "cp %s %s" %(''.join((init_dpff_dir, '/model.*')), model_dir)
    call.call_simple_shell(train_dir, cmd)
  else:
    with open(deepmd_inp_file, 'r') as f:
      deepmd_dic = json.load(f)

    numb_test = deepmd_dic['training']['numb_test']
    init_train_data = deepmd_dic['training']['systems']
    data_dir, final_data_dir = assign_data_dir(work_dir, init_train_data, iter_id, numb_test)

    if ( iter_id > 1 and len(final_data_dir) == 0 ):
      str_print = '''  There are no new data in No.%d iteration, two possible reasons:
  (1) If the accurate ratio is larger than 96.0%%, the deepff is converged, cheers!
  (2) If the accurate ratio is very small (less than 4.0%%), the deepff is not
  converged. But the program have to be stopped. The whole iterations does not
  succeed. The only thing we can do is changing either initial data sets or
  training parameters.''' %(iter_id-1)
      print ('Warning'.center(80,'*'), flush=True)
      print (str_print, flush=True)
      write_data.write_active_data(work_dir, iter_id, tot_atoms_type_dic)
      exit()

    fix_stop_batch = deepmd_test_dic['fix_stop_batch']
    use_prev_model = deepmd_test_dic['use_prev_model']
    lr_scale = deepmd_test_dic['lr_scale']
    start_lr = deepmd_test_dic['start_lr']
    batch_size = deepmd_dic['training']['batch_size']
    if not fix_stop_batch:
      epoch_num = deepmd_test_dic['epoch_num']
      stop_batch = math.ceil(epoch_num*int(sum(data_num)/batch_size[0])/10000)*10000
      decay_steps = math.ceil(int(sum(data_num)/batch_size[0])/1000)*1000
      if ( stop_batch < decay_steps*200 ):
        stop_batch = decay_steps*200
      deepmd_dic['learning_rate']['decay_steps'] = decay_steps
      deepmd_dic['training']['stop_batch'] = stop_batch
    else:
      deepmd_dic['training']['stop_batch'] = deepmd_test_dic['stop_batch']

    deepmd_dic['training']['systems'] = data_dir
    deepmd_dic['training']['batch_size'] = [batch_size[0]]*len(data_dir)

    if use_prev_model:
      if ( lr_scale**iter_id > 100.0 ):
        deepmd_dic['learning_rate']['start_lr'] = start_lr/(100.0)
      else:
        deepmd_dic['learning_rate']['start_lr'] = start_lr/(lr_scale**iter_id)
      deepmd_dic['training']['stop_batch'] = int(deepmd_test_dic['stop_batch']/2)
      prob_sys_1, prob_sys_2 = assign_prob(data_dir, final_data_dir, data_num, iter_id)
      auto_prob_style = data_op.comb_list_2_str(['prob_sys_size', prob_sys_1, prob_sys_2], ';')
      deepmd_dic['training']['auto_prob_style'] = auto_prob_style

    descr_seed = np.random.randint(10000000000)
    deepmd_dic['model']['descriptor']['seed'] = descr_seed
    fit_seed = np.random.randint(10000000000)
    deepmd_dic['model']['fitting_net']['seed'] = fit_seed
    tra_seed = np.random.randint(10000000000)
    deepmd_dic['training']['seed'] = tra_seed
    json_str = json.dumps(deepmd_dic, indent=4)

    with open(''.join((model_dir, '/input.json')), 'w') as json_file:
      json_file.write(json_str)

if __name__ == '__main__':
  from CP2K_kit.deepff import gen_deepmd_task
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
  gen_deepmd_task.gen_deepmd_task(deepmd_dic, work_dir, 0, init_train_data, seed, numb_test)
  #Test run_deepmd function
  deepmd_run.run_deepmd(work_dir, 0)
