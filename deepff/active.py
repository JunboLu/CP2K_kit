#! /usr/env/bin python

import os
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import read_input
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import file_tools
from CP2K_kit.deepff import check_deepff
from CP2K_kit.deepff import load_data
from CP2K_kit.deepff import deepmd_run
from CP2K_kit.deepff import lammps_run
from CP2K_kit.deepff import model_devi
from CP2K_kit.deepff import cp2k_run
from CP2K_kit.deepff import sysinfo

def dump_input(work_dir, inp_file, f_key):

  '''
  dump_input: dump deepff input file, it will call read_input module.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    inp_file: string
      inp_file is the deepff input file
    f_key: 1-d string list
      f_key is fixed to: ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  Returns :
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic: dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic: dictionary
      cp2k_dic contains keywords used in cp2k.
    model_devi_dic: dictionary
      model_devi contains keywords used in model_devi.
    environ_dic: dictionary
      environ_dic contains keywords used in environment.
  '''

  job_type_param = read_input.dump_info(work_dir, inp_file, f_key)
  deepmd_dic = job_type_param[0]
  lammps_dic = job_type_param[1]
  cp2k_dic = job_type_param[2]
  model_devi_dic = job_type_param[3]
  environ_dic = job_type_param[4]

  return deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic

def get_atoms_type(deepmd_dic):

  '''
  get_atoms_type: get atoms type for total systems

  Args:
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
  Returns:
    final_atoms_type: list
      final_atoms_type is the atoms type for all systems.
      Example: ['O', 'H']
  '''

  import linecache

  atoms_type = []
  train_dic = deepmd_dic['training']
  for key in train_dic:
    if ( 'system' in key ):
      traj_coord_file = train_dic[key]['traj_coord_file']
      line_num = file_tools.grep_line_num("'PDB file'", traj_coord_file, os.getcwd())
      if ( line_num == 0 ):
        coord_file_type = 'coord_xyz'
      else:
        coord_file_type = 'coord_pdb'
      atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_id, end_id, time_step = \
      traj_info.get_traj_info(traj_coord_file, coord_file_type)
      atoms = []
      for i in range(atoms_num):
        line_i = linecache.getline(traj_coord_file, pre_base+pre_base_block+i+1)
        line_i_split = data_op.split_str(line_i, ' ', '\n')
        if ( coord_file_type == 'coord_xyz' ):
          atoms.append(line_i_split[0])
        elif ( coord_file_type == 'coord_pdb' ):
          atoms.append(line_i_split[len(line_i_split)-1])
      linecache.clearcache()
      atoms_type.append(data_op.list_replicate(atoms))

  tot_atoms_type = data_op.list_reshape(atoms_type)
  final_atoms_type = data_op.list_replicate(tot_atoms_type)

  return final_atoms_type

def dump_init_data(work_dir, deepmd_dic, train_stress, tot_atoms_type_dic):

  '''
  dump_init_data: load initial training data.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    train_stress: bool
      train_stress is whether we need to dump stress.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    init_train_data: 1-d string list
      init_train_data contains initial training data directories.
    init_data_num : int
      init_data_num is the number of data for initial training.
  '''

  init_train_data_dir = ''.join((work_dir, '/init_train_data'))
  if ( not os.path.exists(init_train_data_dir) ):
    cmd = "mkdir %s" % ('init_train_data')
    call.call_simple_shell(work_dir, cmd)

  i = 0
  init_train_data = []
  init_data_num = 0
  train_dic = deepmd_dic['training']
  shuffle_data = train_dic['shuffle_data']
  for key in train_dic:
    if ( 'system' in key):
      save_dir = ''.join((work_dir, '/init_train_data/data_', str(i)))
      if ( not os.path.exists(save_dir) ):
        cmd = "mkdir %s" % (save_dir)
        call.call_simple_shell(work_dir, cmd)
      init_train_data.append(save_dir)
      cmd = "ls | grep %s" %("'set.'")
      set_dir_name = call.call_returns_shell(save_dir, cmd)
      choosed_num = train_dic[key]['choosed_frame_num']
      data_num = []
      if ( len(set_dir_name) > 0 ):
        for set_dir in set_dir_name:
          data_num_part = []
          set_dir_abs = ''.join((save_dir, '/', set_dir))
          coord_npy_file = ''.join((set_dir_abs, '/coord.npy'))
          force_npy_file = ''.join((set_dir_abs, '/force.npy'))
          box_npy_file = ''.join((set_dir_abs, '/box.npy'))
          energy_npy_file = ''.join((set_dir_abs, '/energy.npy'))
          if ( all(os.path.exists(npy_file) for npy_file in [coord_npy_file, force_npy_file, box_npy_file, energy_npy_file]) ):
            for npy_file in [coord_npy_file, force_npy_file, box_npy_file, energy_npy_file]:
              data_num_part.append(len(np.load(npy_file)))
          else:
            data_num_part = [0,0,0,0]
          virial_npy_file = ''.join((set_dir_abs, '/virial.npy'))
          if ( os.path.exists(virial_npy_file) ):
            data_num_part.append(len(np.load(virial_npy_file)))
          data_num.append(data_num_part)
      else:
        data_num = [[0,0,0,0]]
      data_num = data_op.add_2d_list(data_num)
      if ( all(j == choosed_num for j in data_num) ):
        if ( len(set_dir_name) == 1 ):
          init_data_num_part = choosed_num
        else:
          final_set_dir_abs = ''.join((save_dir, '/', set_dir_name[len(set_dir_name)-1]))
          final_energy_npy_file = ''.join((final_set_dir_abs, '/energy.npy'))
          init_data_num_part = choosed_num-len(np.load(final_energy_npy_file))
      else:
        traj_type = train_dic[key]['traj_type']
        start = train_dic[key]['start_frame']
        end = train_dic[key]['end_frame']
        parts = train_dic[key]['set_parts']
        if ( traj_type == 'md' ):
          traj_coord_file = train_dic[key]['traj_coord_file']
          traj_frc_file = train_dic[key]['traj_frc_file']
          traj_cell_file = train_dic[key]['traj_cell_file']
          traj_stress_file = train_dic[key]['traj_stress_file']
          load_data.load_data_from_dir(traj_coord_file, traj_frc_file, traj_cell_file, traj_stress_file, \
                                       train_stress, work_dir, save_dir, start, end, choosed_num, tot_atoms_type_dic)
        elif ( traj_type == 'mtd' ):
          data_dir = train_dic[key]['data_dir']
          task_dir_prefix = train_dic[key]['task_dir_prefix']
          proj_name = train_dic[key]['proj_name']
          out_file_name = train_dic[key]['out_file_name']
          load_data.load_data_from_sepfile(data_dir, task_dir_prefix, proj_name, tot_atoms_type_dic, \
                                           save_dir, start, end, choosed_num, out_file_name)
        energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(save_dir)
        init_data_num_part, init_test_data_num_part = load_data.raw_data_to_set(parts, shuffle_data, save_dir, energy_array, \
                                                                                coord_array, frc_array, box_array, virial_array)
      init_data_num = init_data_num+init_data_num_part
      i = i+1

  if ( 'set_data_dir' in train_dic.keys() ):
    init_train_data.append(os.path.abspath(train_dic['set_data_dir']))
    energy_npy_file = ''.join((os.path.abspath(train_dic['set_data_dir']), '/set.000/energy.npy'))
    set_data_num = len(np.load(energy_npy_file))
    init_data_num = init_data_num+set_data_num

  return init_train_data, init_data_num

def run_iter(inp_file, deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, \
             init_train_data, init_data_num, work_dir, tot_atoms_type_dic):

  '''
  run_iter: run active learning iterations.

  Args:
    inp_file: string
      inp_file is the input file of CP2K_kit.
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic: dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic: dictionary
      cp2k_dic contains keywords used in cp2k.
    model_devi_dic: dictionary
      model_devi_dic contains keywords used in model_devi.
    environ_dic: dictionary
      environ_dic contains keywords used in environment.
    init_train_data: 1-d string list
      init_train_data contains initial training data directories.
    init_data_num: int
      init_data_num is the data number for initial training.
    work_dir: string
      work_dir is working directory of CP2K_kit.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  proc_num, proc_num_per_node, host, ssh = sysinfo.get_host(work_dir)
  device, usage = sysinfo.analyze_gpu(host, ssh, work_dir)

  max_iter = model_devi_dic['max_iter']
  restart_iter = model_devi_dic['restart_iter']

  if ( restart_iter == 0 ):
    data_num = [init_data_num]
  else:
    data_num = model_devi_dic['data_num']
  restart_stage = model_devi_dic['restart_stage']

  cmd = "ls | grep %s" %("'iter_'")
  iter_info = call.call_returns_shell(work_dir, cmd)
  if ( restart_iter == 0 and restart_stage == 0 and len(iter_info) > 0 ):
    iter_0_dir = ''.join((work_dir, '/iter_0'))
    dir_num = len(call.call_returns_shell(iter_0_dir, "ls -ll |awk '/^d/ {print $NF}'"))
    if ( dir_num > 1 ):
      log_info.log_error('There are iteration directories in %s, please use CP2K_KIT.restart file in %s as input.' \
                          %(work_dir, work_dir), 'Warning')
      exit()

  numb_test = deepmd_dic['training']['numb_test']
  model_type = deepmd_dic['training']['model_type']
  neuron = deepmd_dic['training']['neuron']
  shuffle_data = deepmd_dic['training']['shuffle_data']
  train_stress = deepmd_dic['training']['train_stress']
  use_prev_model = deepmd_dic['training']['use_prev_model']
  lr_scale = deepmd_dic['training']['lr_scale']

  nsteps = int(lammps_dic['nsteps'])
  change_init_str = lammps_dic['change_init_str']

  model_devi_freq = int(model_devi_dic['model_devi_freq'])
  conv_new_data_num = int(nsteps/model_devi_freq*0.04)
  choose_new_data_num_limit = model_devi_dic['choose_new_data_num_limit']
  force_conv = model_devi_dic['force_conv']
  model_devi_steps = int(nsteps/model_devi_freq)+1

  cp2k_exe = environ_dic['cp2k_exe']
  cp2k_env_file = environ_dic['cp2k_env_file']
  parallel_exe = environ_dic['parallel_exe']
  cuda_dir = environ_dic['cuda_dir']
  cp2k_job_per_node = environ_dic['cp2k_job_per_node']
  lmp_job_per_node = environ_dic['lmp_job_per_node']
  lmp_mpi_num_per_job = environ_dic['lmp_mpi_num_per_job']
  lmp_omp_num_per_job = environ_dic['lmp_omp_num_per_job']

  dp_path = sysinfo.get_dp_path(work_dir)
  lmp_exe, lmp_path = sysinfo.get_lmp_path(work_dir)
  mpi_path = sysinfo.get_mpi_path(work_dir)

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
      check_deepff.write_restart_inp(inp_file, i, 0, data_num, work_dir)
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

      deepmd_run.gen_deepmd_task(deepmd_dic, work_dir, i, init_train_data, numb_test, \
                                 descr_seed, fit_seed, tra_seed, neuron, model_type, data_num, lr_scale)
      deepmd_run.run_deepmd(work_dir, i, use_prev_model, parallel_exe, dp_path, host, device, usage, cuda_dir)
      check_deepff.write_restart_inp(inp_file, i, 1, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 ):
      #Perform lammps calculations
      print ('Step 2: lammps tasks', flush=True)

      lammps_run.gen_lmpmd_task(lammps_dic, work_dir, i, change_init_str, tot_atoms_type_dic)
      lammps_run.run_lmpmd(work_dir, i, lmp_path, lmp_exe, parallel_exe, mpi_path, lmp_job_per_node, \
                           lmp_mpi_num_per_job, lmp_omp_num_per_job, proc_num_per_node, host, ssh, device)
      check_deepff.write_restart_inp(inp_file, i, 2, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 ):
      #Perform lammps force calculations
      if ( restart_stage == 2 ):
        print ('Step 2: lammps tasks', flush=True)
      sys_num, atoms_type_multi_sys, atoms_num_tot = lammps_run.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      lammps_run.gen_lmpfrc_file(work_dir, i, atoms_num_tot, atoms_type_multi_sys)
      lammps_run.run_lmpfrc(work_dir, i, lmp_path, lmp_exe, mpi_path, parallel_exe, \
                            proc_num_per_node, host, ssh, atoms_num_tot)
      check_deepff.write_restart_inp(inp_file, i, 3, data_num, work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 or restart_stage == 3 ):
      #Get force-force correlation and then choose new structures
      sys_num, atoms_type_multi_sys, atoms_num_tot = lammps_run.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      struct_index, success_ratio_sys, success_ratio, success_devi_ratio = \
      model_devi.choose_lmp_str(work_dir, i, atoms_type_multi_sys, force_conv)

      for j in range(len(success_ratio_sys)):
        print ('  The accurate ratio for system %d in iteration %d is %.2f%%' %(j, i, success_ratio_sys[j]*100), flush=True)

      print ('  The accurate ratio for whole %d systems in iteration %d is %.2f%%' \
             %(sys_num, i, success_ratio*100), flush=True)
      print ('  The accurate deviation ratio for whole %d systems in iteration %d is %.2f%%' \
             %(sys_num, i, success_devi_ratio*100), flush=True)

      if ( min(success_ratio_sys) >= 0.96 and success_ratio+success_devi_ratio >= 0.999 ):
        print (''.center(80,'*'), flush=True)
        print ('Cheers! deepff is converged!', flush=True)
        if ( i != 0 ):
          check_deepff.write_active_data(work_dir, i, tot_atoms_type_dic)
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

      print ('Step 3: cp2k tasks', flush=True)
      #Perform cp2k calculation
      cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, i, atoms_type_multi_sys, atoms_num_tot, struct_index, \
                             conv_new_data_num, choose_new_data_num_limit, train_stress, success_ratio, success_devi_ratio)
      cp2k_run.run_cp2kfrc(work_dir, i, cp2k_exe, parallel_exe, cp2k_env_file, \
                           cp2k_job_per_node, proc_num_per_node, host, ssh, atoms_num_tot)

      #Dump new data of cp2k
      for j in range(sys_num):
        file_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc/sys_', str(j)))
        load_data.load_data_from_sepfile(file_dir, 'task_', 'cp2k', tot_atoms_type_dic)
        cp2k_data_dir = ''.join((file_dir, '/data'))
        if ( os.path.exists(cp2k_data_dir) ):
          energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(cp2k_data_dir)
          train_data_num, test_data_num = load_data.raw_data_to_set(1, shuffle_data, cp2k_data_dir, energy_array, \
                                                                    coord_array, frc_array, box_array, virial_array)
          if ( test_data_num > numb_test ):
            data_num.append(train_data_num)
          model_devi_steps = int(nsteps/model_devi_freq)+1
          if ( test_data_num < numb_test and success_ratio < float((model_devi_steps-train_data_num)/model_devi_steps) ):
            log_info.log_error('Warning: little selected structures, check the deepmd training.')
            exit()

      print ('  Success: dump new raw data of cp2k', flush=True)
      restart_stage = 0

    check_deepff.write_restart_inp(inp_file, i+1, 0, data_num, work_dir)

    if ( i == max_iter-1 ):
      log_info.log_error('Active learning does not converge')
      check_deepff.write_active_data(work_dir, i+1, tot_atoms_type_dic)

def kernel(work_dir, inp_file):

  '''
  kernel: kernel function to do active learning.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    inp_file: string
      inp_file is the deepff input file
  '''

  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']

  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  dump_input(work_dir, inp_file, deepff_key)

  proc_num, proc_num_per_node, host, ssh = sysinfo.get_host(work_dir)
  if ( len(data_op.list_replicate(proc_num_per_node)) != 1 ):
    host_str = data_op.comb_list_2_str(host, ' ')
    log_info.log_error('Resource error: the number of cores in %s are not equal, please submit your jobs to nodes with same number of cores' %(host_str))
    exit()

  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num_per_node[0])

  train_stress = deepmd_dic['training']['train_stress']

  tot_atoms_type = get_atoms_type(deepmd_dic)
  tot_atoms_type_dic = OrderedDict()
  for i in range(len(tot_atoms_type)):
    tot_atoms_type_dic[tot_atoms_type[i]] = i

  if ( deepmd_dic['model']['type_map'] != tot_atoms_type ):
    type_map_str = data_op.comb_list_2_str(tot_atoms_type, ' ')
    log_info.log_error('Input error: type_map should be %s, please reset deepff/deepmd/model/type_map' %(type_map_str))
    exit()

  if ( len(deepmd_dic['model']['descriptor']['sel']) != len(tot_atoms_type) ):
    log_info.log_error('Input error: sel should be %d integers, please reset deepff/deepmd/model/descriptor/sel' %(len(tot_atoms_type)))
    exit()

  init_train_data, init_data_num = dump_init_data(work_dir, deepmd_dic, train_stress, tot_atoms_type_dic)

  print ('Check input file: no error in %s' %(inp_file), flush=True)
  print ('Initial training data:', flush=True)
  for i in range(len(init_train_data)):
    print ('%s' %(data_op.str_wrap(init_train_data[i], 80)), flush=True)

  run_iter(inp_file, deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, \
           init_train_data, init_data_num, work_dir, tot_atoms_type_dic)

if __name__ == '__main__':

  from CP2K_kit.deepff import active
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  inp_file = 'input.inp'
  max_cycle = 100
  active.kernel(work_dir, inp_file)
