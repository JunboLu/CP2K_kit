#! /usr/env/bin python

import os
import json
import copy
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.deepff import check_deepff
from CP2K_kit.deepff import load_data
from CP2K_kit.deepff import gen_deepmd_task
from CP2K_kit.deepff import deepmd_run
from CP2K_kit.deepff import gen_lammps_task
from CP2K_kit.deepff import lammps_run
from CP2K_kit.deepff import model_devi
from CP2K_kit.deepff import dp_test
from CP2K_kit.deepff import gen_cp2k_task
from CP2K_kit.deepff import cp2k_run
from CP2K_kit.deepff import sys_info
from CP2K_kit.deepff import write_data
from CP2K_kit.deepff import process

def model_devi_iter(work_dir, inp_file, deepmd_dic, lammps_dic, cp2k_dic, active_learn_dic, environ_dic, \
                    init_train_data, init_data_num, tot_atoms_type_dic):

  '''
  model_devi_iter: run active learning iterations with model deviation.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    inp_file: string
      inp_file is the input file of CP2K_kit.
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic: dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic: dictionary
      cp2k_dic contains keywords used in cp2k.
    active_learn_dic: dictionary
      active_learn_dic contains keywords used in active learning.
    environ_dic: dictionary
      environ_dic contains keywords used in environment.
    init_train_data: 1-d string list
      init_train_data contains initial training data directories.
    init_data_num: int
      init_data_num is the data number for initial training.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  proc_num, proc_num_per_node, host, ssh = sys_info.get_host(work_dir)
  device = sys_info.analyze_gpu(host, ssh, work_dir)
  device_num = []
  for i in device:
    device_num.append(len(i))
  if ( len(data_op.list_replicate(device_num)) != 1 and 0 in device_num ):
    log_info.log_error("Resource error: we recommend users use pure CPU or GPU nodes")
    exit()

  max_iter = active_learn_dic['max_iter']
  restart_iter = active_learn_dic['restart_iter']

  if ( restart_iter == 0 ):
    data_num = [init_data_num]
  else:
    data_num = active_learn_dic['data_num']
  restart_stage = active_learn_dic['restart_stage']

  cmd = "ls | grep %s" %("'iter_'")
  iter_info = call.call_returns_shell(work_dir, cmd)
  if ( restart_iter == 0 and restart_stage == 0 and len(iter_info) > 0 ):
    iter_0_dir = ''.join((work_dir, '/iter_0'))
    dir_num = len(call.call_returns_shell(iter_0_dir, "ls -ll |awk '/^d/ {print $NF}'"))
    if ( dir_num > 1 ):
      log_info.log_error('There are iteration directories in %s, please use CP2K_KIT.restart file in %s as input.' \
                          %(work_dir, work_dir), 'Warning')
      exit()

  atom_mass_dic = deepmd_dic['model']['atom_mass']
  numb_test = deepmd_dic['training']['numb_test']
  model_type = deepmd_dic['training']['model_type']
  neuron = deepmd_dic['training']['neuron']
  shuffle_data = deepmd_dic['training']['shuffle_data']
  train_stress = deepmd_dic['training']['train_stress']
  use_prev_model = deepmd_dic['training']['use_prev_model']

  nsteps = int(lammps_dic['nsteps'])
  judge_freq = int(active_learn_dic['judge_freq'])
  conv_new_data_num = int(nsteps/judge_freq*0.04)
  choose_new_data_num_limit = active_learn_dic['choose_new_data_num_limit']
  success_force_conv = active_learn_dic['success_force_conv']
  max_force_conv = active_learn_dic['max_force_conv']
  active_learn_steps = int(nsteps/judge_freq)+1

  cp2k_exe = environ_dic['cp2k_exe']
  cp2k_env_file = environ_dic['cp2k_env_file']
  parallel_exe = environ_dic['parallel_exe']
  cuda_dir = environ_dic['cuda_dir']
  dp_version = environ_dic['dp_version']
  cp2k_job_per_node = environ_dic['cp2k_job_per_node']
  lmp_job_per_node = environ_dic['lmp_job_per_node']
  dp_job_per_node = environ_dic['dp_job_per_node']
  lmp_mpi_num_per_job = environ_dic['lmp_mpi_num_per_job']
  lmp_omp_num_per_job = environ_dic['lmp_omp_num_per_job']

  dp_path = sys_info.get_dp_path(work_dir)
  lmp_exe, lmp_path = sys_info.get_lmp_path(work_dir)
  mpi_path = sys_info.get_mpi_path(work_dir)

  #np.random.seed(1234567890)

  for i in range(restart_iter, max_iter, 1):

    print (''.join(('iter_', str(i))).center(80,'*'), flush=True)

    #Generate iteration directory
    iter_restart = ''.join(('iter_', str(i)))
    iter_restart_dir = ''.join((work_dir, '/', iter_restart))

    if ( not os.path.exists(iter_restart_dir) ):
      cmd = "mkdir %s" % (iter_restart)
      call.call_simple_shell(work_dir, cmd)

    if ( restart_stage == 0 ):
      #Perform deepmd calculation
      write_data.write_restart_inp(inp_file, i, 0, data_num, work_dir)
      print ('Step 1: deepmd-kit tasks', flush=True)

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
          descr_seed.append(np.random.randint(10000000000))
          fit_seed.append(np.random.randint(10000000000))
          tra_seed.append(np.random.randint(10000000000))

      if ( model_type == 'use_node' ):
        descr_seed = []
        fit_seed = []
        tra_seed = []

        for j in range(len(neuron)):
          descr_seed.append(np.random.randint(10000000000))
          fit_seed.append(np.random.randint(10000000000))
          tra_seed.append(np.random.randint(10000000000))

      gen_deepmd_task.gen_deepmd_model_task(deepmd_dic, work_dir, i, init_train_data, numb_test, descr_seed, fit_seed, \
                                            tra_seed, neuron, model_type, data_num, tot_atoms_type_dic)
      deepmd_run.run_deepmd(work_dir, i, use_prev_model, parallel_exe, dp_path, host, ssh, device, cuda_dir, dp_version, dp_job_per_node)
      write_data.write_restart_inp(inp_file, i, 1, data_num, work_dir)
      failure_model = process.check_deepff_run(work_dir, i)
      if ( len(failure_model) == 0 ):
        pass
      else:
        print ('Warning'.center(80,'*'), flush=True)
        for model_id in failure_model:
          str_print = 'The training is failed as force is fluctuated in No.%d model in No.%d iteration' %(model_id, i)
          str_print = data_op.str_wrap(str_print, 80, '')
          print (str_print, flush=True)
        exit()

    if ( restart_stage == 0 or restart_stage == 1 ):
      #Perform lammps calculations
      print ('Step 2: lammps tasks', flush=True)

      gen_lammps_task.gen_lmpmd_task(lammps_dic, work_dir, i, atom_mass_dic, tot_atoms_type_dic)
      lammps_run.run_lmpmd(work_dir, i, lmp_path, lmp_exe, parallel_exe, mpi_path, lmp_job_per_node, \
                           lmp_mpi_num_per_job, lmp_omp_num_per_job, proc_num_per_node, host, ssh, device)
      write_data.write_restart_inp(inp_file, i, 2, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 ):
      #Perform lammps force calculations
      if ( restart_stage == 2 ):
        print ('Step 2: lammps tasks', flush=True)
      sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot = process.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      gen_lammps_task.gen_lmpfrc_file(work_dir, i, atom_mass_dic, atoms_num_tot, atoms_type_multi_sys, use_mtd_tot, 'model_devi')
      lammps_run.run_lmpfrc(work_dir, i, lmp_path, lmp_exe, mpi_path, parallel_exe, \
                            proc_num_per_node, host, ssh, atoms_num_tot, use_mtd_tot, 'model_devi')
      write_data.write_restart_inp(inp_file, i, 3, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 or restart_stage == 3 ):
      #Get force-force correlation and then choose new structures
      print ('step 3: model deviation', flush=True)
      sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot = process.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      struct_index, success_ratio_sys, success_ratio, success_devi_ratio = \
      model_devi.choose_lmp_str(work_dir, i, atoms_type_multi_sys, success_force_conv, max_force_conv)

      for j in range(len(success_ratio_sys)):
        print ('  The accurate ratio for system %d in iteration %d is %.2f%%' %(j, i, success_ratio_sys[j]*100), flush=True)

      print ('  The accurate ratio for whole %d systems in iteration %d is %.2f%%' \
             %(sys_num, i, success_ratio*100), flush=True)
      print ('  The accurate deviation ratio for whole %d systems in iteration %d is %.2f%%' \
             %(sys_num, i, success_devi_ratio*100), flush=True)

      if ( min(success_ratio_sys) >= 0.95 and success_ratio+success_devi_ratio > 0.99 ):
        print (''.center(80,'*'), flush=True)
        print ('Cheers! deepff is converged!', flush=True)
        if ( i != 0 ):
          write_data.write_active_data(work_dir, i, tot_atoms_type_dic)
        exit()

      total_task_num = 0
      for sys_id in struct_index:
        total_task_num_sys = 0
        for task_id in struct_index[sys_id]:
          total_task_num_sys = total_task_num_sys+len(struct_index[sys_id][task_id])
        total_task_num = total_task_num+total_task_num_sys
      if ( total_task_num == 0 ):
        log_info.log_error('Warning: No selected structure for cp2k calculations, check the deepmd training.')
        exit()

      print ('Step 4: cp2k tasks', flush=True)
      #Perform cp2k calculation
      gen_cp2k_task.gen_cp2k_task(cp2k_dic, work_dir, i, atoms_type_multi_sys, atoms_num_tot, struct_index, conv_new_data_num, \
                                  choose_new_data_num_limit, train_stress, 'model_devi', success_ratio, success_devi_ratio)
      cp2k_run.run_cp2kfrc(work_dir, i, cp2k_exe, parallel_exe, cp2k_env_file, \
                           cp2k_job_per_node, proc_num_per_node, host, ssh, atoms_num_tot)

      #Dump new data of cp2k
      for j in range(sys_num):
        cp2k_sys_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc/sys_', str(j)))
        task_num, task_dir = process.get_task_num(cp2k_sys_dir, True)
        for k in range(task_num):
          cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/', task_dir[k]))
          traj_num = process.get_traj_num(cp2k_sys_task_dir)
          if ( traj_num != 0 ):
            data_dir = ''.join((cp2k_sys_task_dir, '/data'))
            if ( not os.path.exists(data_dir) ):
              cmd = "mkdir %s" % ('data')
              call.call_simple_shell(cp2k_sys_task_dir, cmd)
            total_index = data_op.gen_list(0, traj_num-1, 1)
            total_index_array = np.array(total_index)
            np.random.shuffle(total_index_array)
            choosed_index = list(total_index_array[0:traj_num])
            load_data.load_data_from_sepfile(cp2k_sys_task_dir, data_dir, 'traj_', 'cp2k', tot_atoms_type_dic, choosed_index)
            energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(data_dir)
            train_data_num, test_data_num = load_data.raw_data_to_set(1, shuffle_data, data_dir, energy_array, \
                                                                    coord_array, frc_array, box_array, virial_array)
            if ( test_data_num > numb_test ):
              data_num.append(train_data_num)
            if ( test_data_num < numb_test and success_ratio < float((active_learn_steps-train_data_num)/active_learn_steps) ):
              log_info.log_error('Warning: little selected structures, check the deepmd training.')
              exit()

      print ('  Success: dump new raw data of cp2k', flush=True)
      restart_stage = 0

    write_data.write_restart_inp(inp_file, i+1, 0, data_num, work_dir)

    if ( i == max_iter-1 ):
      log_info.log_error('Active learning does not converge')
      write_data.write_active_data(work_dir, i+1, tot_atoms_type_dic)

def dp_test_iter(work_dir, inp_file, deepmd_dic, lammps_dic, active_learn_dic, cp2k_dic, environ_dic):

  '''
  dp_test_iter: run active learning iterations with deepmd test.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    inp_file: string
      inp_file is the input file of CP2K_kit.
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic: dictionary
      lammpd_dic contains keywords used in lammps.
    active_learn_dic: dictionary
      active_learn_dic contains keywords used in active learning.
    cp2k_dic: dictionary
      cp2k_dic contains keywords used in cp2k.
    environ_dic: dictionary
      environ_dic contains keywords used in environment.
  Returns:
    none
  '''

  proc_num, proc_num_per_node, host, ssh = sys_info.get_host(work_dir)
  device = sys_info.analyze_gpu(host, ssh, work_dir)
  device_num = []
  for i in device:
    device_num.append(len(i))
  if ( len(data_op.list_replicate(device_num)) != 1 and 0 in device_num ):
    log_info.log_error("Resource error: we recommend users use pure CPU or GPU nodes")
    exit()

  max_iter = active_learn_dic['max_iter']
  restart_iter = active_learn_dic['restart_iter']
  data_num = active_learn_dic['data_num']
  restart_stage = active_learn_dic['restart_stage']

  cmd = "ls | grep %s" %("'iter_'")
  iter_info = call.call_returns_shell(work_dir, cmd)
  if ( restart_iter == 0 and restart_stage == 0 and len(iter_info) > 0 ):
    iter_0_dir = ''.join((work_dir, '/iter_0'))
    dir_num = len(call.call_returns_shell(iter_0_dir, "ls -ll |awk '/^d/ {print $NF}'"))
    if ( dir_num > 1 ):
      log_info.log_error('There are iteration directories in %s, please use CP2K_KIT.restart file in %s as input.' \
                          %(work_dir, work_dir), 'Warning')
      exit()

  init_dpff_dir = deepmd_dic['init_dpff_dir']
  use_prev_model = deepmd_dic['use_prev_model']
  train_stress = deepmd_dic['train_stress']
  deepmd_inp_file = ''.join((init_dpff_dir, '/input.json'))
  if ( not os.path.exists(deepmd_inp_file) ):
    log_info.log_error("Running error: %s does not exist" %(deepmd_inp_file))
    exit()

  with open(deepmd_inp_file, 'r') as f:
    deepmd_dic_json = json.load(f)

  atom_mass_dic = deepmd_dic['atom_mass']
  shuffle_data = deepmd_dic['shuffle_data']
  numb_test = deepmd_dic_json['training']['numb_test']
  tot_atoms_type = deepmd_dic_json['model']['type_map']
  tot_atoms_type_dic = OrderedDict()
  for i in range(len(tot_atoms_type)):
    tot_atoms_type_dic[tot_atoms_type[i]] = i

  nsteps = int(lammps_dic['nsteps'])
  change_init_str = lammps_dic['change_init_str']
  judge_freq = int(active_learn_dic['judge_freq'])
  conv_new_data_num = int(nsteps/judge_freq*0.04)
  choose_new_data_num_limit = active_learn_dic['choose_new_data_num_limit']
  data_num = active_learn_dic['data_num']
  success_force_conv = active_learn_dic['success_force_conv']
  max_force_conv = active_learn_dic['max_force_conv']
  energy_conv = active_learn_dic['energy_conv']

  cp2k_exe = environ_dic['cp2k_exe']
  cp2k_env_file = environ_dic['cp2k_env_file']
  parallel_exe = environ_dic['parallel_exe']
  cuda_dir = environ_dic['cuda_dir']
  dp_version = environ_dic['dp_version']
  cp2k_job_per_node = environ_dic['cp2k_job_per_node']
  lmp_job_per_node = environ_dic['lmp_job_per_node']
  dp_job_per_node = environ_dic['dp_job_per_node']
  lmp_mpi_num_per_job = environ_dic['lmp_mpi_num_per_job']
  lmp_omp_num_per_job = environ_dic['lmp_omp_num_per_job']

  dp_path = sys_info.get_dp_path(work_dir)
  lmp_exe, lmp_path = sys_info.get_lmp_path(work_dir)
  mpi_path = sys_info.get_mpi_path(work_dir)

  for i in range(restart_iter, max_iter, 1):

    print (''.join(('iter_', str(i))).center(80,'*'), flush=True)

    #Generate iteration directory
    iter_restart = ''.join(('iter_', str(i)))
    iter_restart_dir = ''.join((work_dir, '/', iter_restart))

    if ( not os.path.exists(iter_restart_dir) ):
      cmd = "mkdir %s" % (iter_restart)
      call.call_simple_shell(work_dir, cmd)

    if ( restart_stage == 0 ):
      #Perform deepmd calculation
      write_data.write_restart_inp(inp_file, i, 0, data_num, work_dir)
      print ('Step 1: deepmd-kit tasks', flush=True)

      gen_deepmd_task.gen_deepmd_test_task(deepmd_dic, work_dir, i, data_num, tot_atoms_type_dic)
      if ( i>0 ):
        deepmd_run.run_deepmd(work_dir, i, use_prev_model, parallel_exe, dp_path, host, ssh, device, cuda_dir, dp_version, dp_job_per_node)
      else:
        str_print = 'Success: the initial deep potential file is copied in %s' %(''.join((work_dir, '/iter_0/01.train/0')))
        str_print = data_op.str_wrap(str_print, 80, '  ')
        print (str_print, flush=True)
      write_data.write_restart_inp(inp_file, i, 1, data_num, work_dir)
      if ( i>0 ):
        failure_model = process.check_deepff_run(work_dir, i)
        if ( len(failure_model) == 0 ):
          pass
        else:
          print ('Warning'.center(80,'*'), flush=True)
          for model_id in failure_model:
            str_print = 'The training is failed as force is fluctuated in No.%d model in No.%d iteration' %(model_id, i)
            str_print = data_op.str_wrap(str_print, 80, '  ')
            print (str_print, flush=True)
          exit()

    if ( restart_stage == 0 or restart_stage == 1 ):
      #Perform lammps calculations
      print ('Step 2: lammps tasks', flush=True)
      gen_lammps_task.gen_lmpmd_task(lammps_dic, work_dir, i, atom_mass_dic, tot_atoms_type_dic)
      lammps_run.run_lmpmd(work_dir, i, lmp_path, lmp_exe, parallel_exe, mpi_path, lmp_job_per_node, \
                           lmp_mpi_num_per_job, lmp_omp_num_per_job, proc_num_per_node, host, ssh, device)
      write_data.write_restart_inp(inp_file, i, 2, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 ):
      #Perform lammps force calculations
      if ( restart_stage == 2 ):
        print ('Step 2: lammps tasks', flush=True)

      sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot = process.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      gen_lammps_task.gen_lmpfrc_file(work_dir, i, atom_mass_dic, atoms_num_tot, atoms_type_multi_sys, use_mtd_tot, 'dp_test')
      lammps_run.run_lmpfrc(work_dir, i, lmp_path, lmp_exe, mpi_path, parallel_exe, \
                            proc_num_per_node, host, ssh, atoms_num_tot, use_mtd_tot, 'dp_test')
      write_data.write_restart_inp(inp_file, i, 3, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 or restart_stage == 3 ):
      print ('Step 3: cp2k tasks', flush=True)
      #Perform cp2k calculation
      sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot = process.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      total_index = OrderedDict()
      for j in range(sys_num):
        total_index_j = OrderedDict()
        lmp_sys_dir = ''.join((work_dir, '/iter_', str(i), '/02.lammps_calc/sys_', str(j)))
        task_num = process.get_task_num(lmp_sys_dir)
        for k in range(task_num):
          lmp_sys_task_data_dir = ''.join((lmp_sys_dir, '/task_', str(k), '/data'))
          data_file_num = process.get_data_num(lmp_sys_task_data_dir)
          total_index_j[k] = range(0, data_file_num, 1)
        total_index[j] = total_index_j
      gen_cp2k_task.gen_cp2k_task(cp2k_dic, work_dir, i, atoms_type_multi_sys, atoms_num_tot, total_index, \
                                  conv_new_data_num, choose_new_data_num_limit, train_stress, 'dp_test')
      cp2k_run.run_cp2kfrc(work_dir, i, cp2k_exe, parallel_exe, cp2k_env_file, \
                           cp2k_job_per_node, proc_num_per_node, host, ssh, atoms_num_tot)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 or restart_stage == 3 or restart_stage == 4 ):
      print ('Step 4: deep potential test', flush=True)
      sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot = process.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      struct_index, success_ratio_sys, success_ratio  = \
      dp_test.active_learning_test(work_dir, i, atoms_type_multi_sys, use_mtd_tot, success_force_conv, max_force_conv, energy_conv)

      for j in range(len(success_ratio_sys)):
        print ('  The accurate ratio for system %d in iteration %d is %.2f%%' %(j, i, success_ratio_sys[j]*100), flush=True)

      print ('  The accurate ratio for whole %d systems in iteration %d is %.2f%%' \
             %(sys_num, i, success_ratio*100), flush=True)

      if ( min(success_ratio_sys) >= 0.95 ):
        print (''.center(80,'*'), flush=True)
        print ('Cheers! deepff is converged!', flush=True)
        if ( i != 0 ):
          write_data.write_active_data(work_dir, i, tot_atoms_type_dic)
        exit()

      #Dump new data of cp2k
      for key in struct_index:
        task_num = len(struct_index[key])

        choosed_task = []
        choosed_index_num = []
        for j in range(task_num):
          choosed_index = struct_index[key][j]
          if ( len(choosed_index) < conv_new_data_num ):
            pass
          else:
            choosed_index_num.append(len(choosed_index))
            choosed_task.append(j)

        choosed_index_num_copy = copy.deepcopy(choosed_index_num)
        if ( sum(choosed_index_num)<choose_new_data_num_limit ):
          pass
        else:
          for j in range(len(choosed_index_num)):
            choosed_index_num[j]=int(choosed_index_num_copy[j]/sum(choosed_index_num_copy)*choose_new_data_num_limit)

        cp2k_sys_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc/sys_', str(key)))
        for j in range(len(choosed_task)):
          choosed_index = struct_index[key][choosed_task[j]]
          choosed_index_array = np.array(choosed_index)
          np.random.shuffle(choosed_index_array)
          choosed_index = list(choosed_index_array[0:choosed_index_num[j]])
          sys_task_index = data_op.comb_list_2_str(sorted(choosed_index), ' ')
          str_print = 'Choosed index for system %d cp2k task %d: %s' %(key, choosed_task[j], sys_task_index)
          str_print = data_op.str_wrap(str_print, 80, '  ')
          print (str_print, flush=True)

          cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(choosed_task[j])))
          data_dir = ''.join((cp2k_sys_task_dir, '/data'))
          if ( not os.path.exists(data_dir) ):
            cmd = "mkdir %s" % ('data')
            call.call_simple_shell(cp2k_sys_task_dir, cmd)
          load_data.load_data_from_sepfile(cp2k_sys_task_dir, data_dir, 'traj_', 'cp2k', tot_atoms_type_dic, sorted(choosed_index))
          energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(data_dir)
          train_data_num, test_data_num = load_data.raw_data_to_set(1, shuffle_data, data_dir, energy_array, \
                                                                    coord_array, frc_array, box_array, virial_array)
          if ( test_data_num > numb_test ):
            data_num.append(train_data_num)
          if ( test_data_num < numb_test and success_ratio < float((active_learn_steps-train_data_num)/active_learn_steps) ):
            log_info.log_error('Warning: little selected structures, check the deepmd training.')
            exit()

      print ('  Success: dump new raw data of cp2k', flush=True)
      restart_stage = 0

    write_data.write_restart_inp(inp_file, i+1, 0, data_num, work_dir)

    if ( i == max_iter-1 ):
      log_info.log_error('Active learning does not converge')
      write_data.write_active_data(work_dir, i+1, tot_atoms_type_dic)

def kernel(work_dir, inp_file, deepff_type):

  '''
  kernel: kernel function to do active learning.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    inp_file: string
      inp_file is the deepff input file.
    deepff_type: string
      deepff_type is the type of deepff.
  Return:
    none
  '''

  proc_num, proc_num_per_node, host, ssh = sys_info.get_host(work_dir)
  if ( len(data_op.list_replicate(proc_num_per_node)) != 1 ):
    host_str = data_op.comb_list_2_str(host, ' ')
    log_info.log_error('Resource error: the number of cores in %s are not equal, please submit your jobs to nodes with same number of cores' %(host_str))
    exit()

  if ( deepff_type == 'active_model_devi' ):
    deepff_key = ['deepmd_model', 'lammps', 'cp2k', 'active_learn', 'environ']
  elif ( deepff_type == 'active_dp_test' ):
    deepff_key = ['deepmd_test', 'lammps', 'cp2k', 'active_learn', 'environ']

  deepmd_dic, lammps_dic, cp2k_dic, active_learn_dic, environ_dic = process.dump_input(work_dir, inp_file, deepff_key)

  if ( deepff_type == 'active_model_devi' ):
    deepmd_dic = check_deepff.check_deepmd_model(deepmd_dic)
  elif ( deepff_type == 'active_dp_test' ):
    deepmd_dic = check_deepff.check_deepmd_test(deepmd_dic)
  active_learn_dic = check_deepff.check_active_learn(active_learn_dic)
  lammps_dic = check_deepff.check_lammps(lammps_dic, active_learn_dic)
  cp2k_dic = check_deepff.check_cp2k(cp2k_dic)
  environ_dic = check_deepff.check_environ(environ_dic, proc_num_per_node[0])

  dp_path = sys_info.get_dp_path(work_dir)
  lmp_exe, lmp_path = sys_info.get_lmp_path(work_dir)
  mpi_path = sys_info.get_mpi_path(work_dir)
  cp2k_exe = environ_dic['cp2k_exe']

  print (data_op.str_wrap('DEEPFF| DEEPMD-KIT EXECUTABLE FILE IS %s' %(dp_path+'/dp'), 80), flush=True)
  print (data_op.str_wrap('DEEPFF| LAMMPS EXECUTABLE FILE IS %s' %(lmp_exe), 80), flush=True)
  if ( deepff_type == 'active_model_devi' ):
    print (data_op.str_wrap('DEEPFF| CP2K EXECUTABLE FILE IS %s' %(cp2k_exe), 80), flush=True)
  elif ( deepff_type == 'active_dp_test' ):
    print (data_op.str_wrap('DEEPFF| CP2K EXECUTABLE FILE IS %s\n' %(cp2k_exe), 80), flush=True)

  if ( deepff_type == 'active_model_devi' ):
    if ( 'set_data_dir' not in deepmd_dic['training'].keys() ):
      tot_atoms_type = process.get_atoms_type(deepmd_dic)
      if ( data_op.reorder_atom_list(deepmd_dic['model']['type_map']) != data_op.reorder_atom_list(tot_atoms_type) ):
        type_map_str = data_op.comb_list_2_str(tot_atoms_type, ' ')
        log_info.log_error('Input error: missing elements in type_map, please reset deepff/deepmd/model/type_map')
        exit()
    else:
      tot_atoms_type = deepmd_dic['model']['type_map']
    tot_atoms_type_dic = OrderedDict()
    for i in range(len(tot_atoms_type)):
      tot_atoms_type_dic[tot_atoms_type[i]] = i

    if ( len(deepmd_dic['model']['descriptor']['sel']) != len(tot_atoms_type) ):
      log_info.log_error('Input error: sel should be %d integers, please reset deepff/deepmd/model/descriptor/sel' %(len(tot_atoms_type)))
      exit()
    train_stress = deepmd_dic['training']['train_stress']
    init_train_data, init_data_num = process.dump_init_data(work_dir, deepmd_dic, train_stress, tot_atoms_type_dic)
    print ('DEEPFF| INITIAL TRAINING DATA:', flush=True)
    for i in range(len(init_train_data)):
      if ( i == len(init_train_data)-1 ):
        print ('%s\n' %(data_op.str_wrap(init_train_data[i], 80)), flush=True)
      else:
        print ('%s' %(data_op.str_wrap(init_train_data[i], 80)), flush=True)

    model_devi_iter(work_dir, inp_file, deepmd_dic, lammps_dic, cp2k_dic, active_learn_dic, \
                    environ_dic, init_train_data, init_data_num, tot_atoms_type_dic)
  elif ( deepff_type == 'active_dp_test' ):
    dp_test_iter(work_dir, inp_file, deepmd_dic, lammps_dic, active_learn_dic, cp2k_dic, environ_dic)

if __name__ == '__main__':

  from CP2K_kit.deepff import active
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  inp_file = 'input.inp'
  max_cycle = 100
  active.kernel(work_dir, inp_file)
