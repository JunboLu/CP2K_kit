#! /usr/env/bin python

import subprocess
import numpy as np
from CP2K_kit.tools import call
from CP2K_kit.tools import read_input
from CP2K_kit.tools import list_dic_op
from CP2K_kit.deepff import load_data
from CP2K_kit.deepff import deepmd_run
from CP2K_kit.deepff import lammps_run
from CP2K_kit.deepff import force_eval
from CP2K_kit.deepff import cp2k_run

def dump_input(work_dir, inp_file, f_key):

  '''
  dump_input : dump deepff input file, it will call read_input module.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    inp_file : string
      inp_file is the deepff input file
    f_key : 1-d string list
      f_key is fixed to: ['deepmd', 'lammps', 'cp2k', 'force_eval', 'environ']
  Returns :
    deepmd_dic : dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic : dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic : dictionary
      cp2k_dic contains keywords used in cp2k.
    force_eval : dictionary
      force_eval contains keywords used in force_eval.
    environ_dic : dictionary
      environ_dic contains keywords used in environment.
  '''

  job_type_param = read_input.dump_info(work_dir, inp_file, f_key)
  deepmd_dic = job_type_param[0]
  lammps_dic = job_type_param[1]
  cp2k_dic = job_type_param[2]
  force_eval_dic = job_type_param[3]
  environ_dic = job_type_param[4]

  return deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic

def dump_init_data(work_dir, deepmd_dic, restart_iter):

  '''
  dump_init_data : load initial training data.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    deepmd_dic : dictionary
      deepmd_dic contains keywords used in deepmd.
    restart_iter : int
      restart_iter is the iteration number of restart.
  Returns :
    init_train_data : 1-d string list
      init_train_data contains initial training data directories.
  '''

  i = 0
  init_train_data = []
  cmd = "mkdir %s" % ('init_train_data')
  if ( restart_iter == 0 ):
    call.call_simple_shell(work_dir, cmd)
  train_dic = deepmd_dic['training']
  for key in train_dic:
    if ( 'system' in key):
      init_train_key_dir = train_dic[key]['directory']
      proj_name = train_dic[key]['proj_name']
      save_dir = ''.join((work_dir, '/init_train_data/data_', str(i)))
      init_train_data.append(save_dir)
      if ( restart_iter == 0 ):
        load_data.load_data_from_dir(init_train_key_dir, work_dir, save_dir, proj_name)
        load_data.raw_to_set(save_dir, 1)
      i = i+1

  return init_train_data

def run_iteration(deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic, \
                  init_train_data, work_dir, max_iter, restart_iter, proc_num, host):

  '''
  run_iteration : run active learning iterations.

  Args :
    deepmd_dic : dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic : dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic : dictionary
      cp2k_dic contains keywords used in cp2k.
    force_eval_dic : dictionary
      force_eval_dic contains keywords used in force_eval.
    environ_dic : dictionary
      environ_dic contains keywords used in environment.
    init_train_data : 1-d string list
      init_train_data contains initial training data directories.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    max_iter : int
      max_iter is the maxium iterations for active learning.
    restart_iter : int
      restart_iter is the iteration number of restart.
    proc_num : int
      proc_num is the number of processor.
    host : 1-d string list
      host is the name of host of computational nodes.
  Returns :
    none
  '''

  import os

  iter_restart = ''.join(('iter_', str(restart_iter)))
  iter_restart_dir = ''.join((work_dir, '/iter_', str(restart_iter)))
  if ( os.path.exists(iter_restart_dir) ):
    cmd = "rm -rf %s" % (iter_restart)
    call.call_simple_shell(work_dir, cmd)

  np.random.seed(12345678)

  for i in range(restart_iter, max_iter, 1):
    #Generate iteration directory
    cmd = "mkdir %s" % (''.join(('iter_', str(i))))
    call.call_simple_shell(work_dir, cmd)

    ##Perform deepmd calculation
    #Get general keywords for deepmd run
    if ( 'numb_test' in deepmd_dic['training'].keys() ):
      numb_test = int(deepmd_dic['training']['numb_test'])
    else:
      numb_test = 5

    if ( 'neuron' in deepmd_dic['training'].keys() ):
      neuron_str = deepmd_dic['training']['neuron']
    else:
      print ('Cannot find neuron, please set neuron in training')
      exit()

    if ( 'model_type' in deepmd_dic['training'].keys() ):
      model_type = deepmd_dic['training']['model_type']
    else:
      model_type = 'use_seed'

    #For different model_type, seed and neuron are different.
    if ( model_type == 'use_seed' ):
      if ( 'seed_num' in deepmd_dic['training'].keys() ):
        seed_num = int(deepmd_dic['training']['seed_num'])
      else:
        seed_num = 4
      descr_seed = []
      fit_seed = []
      tra_seed = []
      for j in range(seed_num):
        descr_seed.append(np.random.randint(100000000))
        fit_seed.append(np.random.randint(100000000))
        tra_seed.append(np.random.randint(100000000))

      neuron = [int(x) for x in neuron_str]

    if ( model_type == 'use_nodes' ):
      neuron = []
      descr_seed = []
      fit_seed = []
      tra_seed = []

      temp_str = ''
      for j in range(len(neuron_str)):
        temp_str = ' '.join((temp_str, neuron_str[j]))
      temp_list = list_dic_op.str_split(temp_str, '...')

      for j in range(len(temp_list)):
        neuron_j = list_dic_op.str_split(temp_list[j], ' ')
        neuron.append([int(x) for x in neuron_j])
        descr_seed.append(np.random.randint(100000000))
        fit_seed.append(np.random.randint(100000000))
        tra_seed.append(np.random.randint(100000000))

    if ( 'parallel_exe' in environ_dic.keys() ):
      parallel_exe = environ_dic['parallel_exe']
    else:
      print ('please set parallel_exe in environ')
      exit()

    if ( 'cuda_dir' in environ_dic.keys() ):
      cuda_dir = environ_dic['cuda_dir']
    else:
      cuda_dir = 'none'

    deepmd_run.gen_deepmd_task(deepmd_dic, work_dir, i, init_train_data, numb_test, \
                               descr_seed, fit_seed, tra_seed, neuron, model_type)
    deepmd_run.run_deepmd(work_dir, i, parallel_exe, host, cuda_dir)

    ##Perform lammps calculation
    atoms_type_dic_tot, atoms_num_tot = \
    lammps_run.gen_lmpmd_task(lammps_dic, work_dir, i)
    sys_num = len(atoms_type_dic_tot)
    lammps_run.run_lmpmd(work_dir, i, proc_num)
    lammps_run.gen_lmpfrc_file(work_dir, i, atoms_num_tot, atoms_type_dic_tot)
    lammps_run.run_lmpfrc(work_dir, i, parallel_exe, proc_num)

    ##Get force-force correlation and then choose new structures
    force_conv = float(force_eval_dic['force_conv'])
    struct_index = force_eval.choose_lmp_str(work_dir, i, atoms_type_dic_tot, atoms_num_tot, force_conv)

    conv_new_data_num = int(force_eval_dic['conv_new_data_num'])
    choose_data_num = []
    for key1 in struct_index:
      for key2 in struct_index[key1]:
        choose_data_num.append(len(struct_index[key1][key2]))
    max_choose_data_num = max(choose_data_num)
    if ( max_choose_data_num <= conv_new_data_num ):
      print ('Cheers! deepff is converged!')
      exit()

    ##Perform cp2k calculation
    if ( 'md_type' in lammps_dic.keys() ):
      md_type = lammps_dic['md_type']
    else:
      md_type = 'nvt'
    if ( md_type == 'nvt' or md_type == 'nve' ):
      get_stress = False
    elif ( md_type == 'npt' ):
      get_stress = True
    choose_new_data_num_limit = int(force_eval_dic['choose_new_data_num_limit'])
    cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, i, atoms_type_dic_tot, atoms_num_tot, \
                           struct_index, conv_new_data_num, choose_new_data_num_limit, get_stress)
    cp2k_run.run_cp2kfrc(work_dir, i, environ_dic, proc_num)

    ##Get new data of CP2K
    for j in range(sys_num):
      file_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc/sys_', str(j)))
      load_data.load_data_from_sepfile(file_dir, 'task_', 'cp2k')
      cp2k_data_dir = ''.join((file_dir, '/data'))
      load_data.raw_to_set(cp2k_data_dir, 1)

def kernel(work_dir, inp_file):

  '''
  kernel : kernel function to do active learning.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    inp_file : string
      inp_file is the deepff input file
  '''

  import os
  import linecache
  import platform
  import multiprocessing

  deepff_key = ['deepmd', 'lammps', 'cp2k', 'force_eval', 'environ']

  deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic = \
  dump_input(work_dir, inp_file, deepff_key)

  if ( 'max_iter' in force_eval_dic.keys() ):
    max_iter = int(force_eval_dic['max_iter'])
  else:
    max_iter = 100

  if ( 'restart_iter' in force_eval_dic.keys() ):
    restart_iter = int(force_eval_dic['restart_iter'])
  else:
    restart_iter = 0

  host_file = ''.join((work_dir, '/hostname'))
  if ( os.path.exists(host_file) ):
    line_num = len(open(host_file).readlines())
    line_1 = linecache.getline(host_file, 1)
    proc_num = int(line_1.strip('\n'))
    host = []
    for i in range(line_num-1):
      line_i = linecache.getline(host_file, i+2)
      line_i_split = list_dic_op.str_split(line_i, ' ')
      host.append(line_i_split[1].strip('\n'))
  else:
    proc_num = int(multiprocessing.cpu_count()/2)
    host = [platform.node()]

  init_train_data = dump_init_data(work_dir, deepmd_dic, restart_iter)

  run_iteration(deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic, \
                init_train_data, work_dir, max_iter, restart_iter, proc_num, host)

if __name__ == '__main__':

  from CP2K_kit.deepff import active
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  inp_file = 'input.inp'
  max_cycle = 100
  active.kernel(work_dir, inp_file)
