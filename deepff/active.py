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

def dump_init_data(work_dir, deepmd_dic, restart_iter, train_stress, tot_atoms_type_dic):

  '''
  dump_init_data: load initial training data.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    deepmd_dic: dictionary
      deepmd_dic contains keywords used in deepmd.
    restart_iter: int
      restart_iter is the iteration number of restart.
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

  if ( restart_iter == 0 ):
    cmd = 'rm -rf init_train_data'
    call.call_simple_shell(work_dir, cmd)

  i = 0
  init_train_data = []
  init_data_num = 0
  cmd = "mkdir %s" % ('init_train_data')
  if ( restart_iter == 0 ):
    call.call_simple_shell(work_dir, cmd)
  train_dic = deepmd_dic['training']
  shuffle_data = train_dic['shuffle_data']
  for key in train_dic:
    if ( 'system' in key):
      save_dir = ''.join((work_dir, '/init_train_data/data_', str(i)))
      init_train_data.append(save_dir)
      if ( restart_iter == 0 ):
        traj_coord_file = train_dic[key]['traj_coord_file']
        traj_frc_file = train_dic[key]['traj_frc_file']
        traj_cell_file = train_dic[key]['traj_cell_file']
        traj_stress_file = train_dic[key]['traj_stress_file']
        start = train_dic[key]['start_frame']
        end = train_dic[key]['end_frame']
        choosed_num = train_dic[key]['choosed_frame_num']
        parts = train_dic[key]['set_parts']
        load_data.load_data_from_dir(traj_coord_file, traj_frc_file, traj_cell_file, traj_stress_file, \
                                     train_stress, work_dir, save_dir, start, end, choosed_num, tot_atoms_type_dic)
        energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(save_dir)
        data_num = load_data.raw_data_to_set(parts, shuffle_data, save_dir, energy_array, coord_array, frc_array, box_array, virial_array)
        init_data_num = init_data_num+data_num
      else:
        if ( os.path.exists(save_dir) ):
          init_data_num = 0
        else:
          log_info.log_error('%s does not exist for training system %d' %(save_dir, i))
          exit()
      i = i+1

  if ( 'set_data_dir' in train_dic.keys() ):
    init_train_data.append(os.path.abspath(train_dic['set_data_dir']))

  return init_train_data, init_data_num

def write_active_data(work_dir, conv_iter, tot_atoms_type_dic):

  '''
  write_active_data: write the data generated by active learning.

  Args:
    work_dir: string
      work_dir is the working directory.
    conv_iter: int
      conv_iter is the number of iteration.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns :
    none
  '''

  active_data_dir = ''.join((work_dir, '/active_data'))

  cmd = "mkdir %s" %('active_data')
  call.call_simple_shell(work_dir, cmd)
  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(''.join((work_dir, '/iter_0/02.lammps_calc')), cmd))

  for i in range(sys_num):

    energy_cp2k = []
    frc_cp2k = []
    frc_x_cp2k = []
    frc_y_cp2k = []
    frc_z_cp2k = []

    sys_dir = ''.join((active_data_dir, '/sys_', str(i)))
    cmd = "mkdir %s" %(''.join(('sys_', str(i))))
    call.call_simple_shell(active_data_dir, cmd)
    energy_file_name = ''.join((sys_dir, '/energy.raw'))
    coord_file_name = ''.join((sys_dir, '/coord.raw'))
    frc_file_name = ''.join((sys_dir, '/force.raw'))
    cell_file_name = ''.join((sys_dir, '/box.raw'))
    energy_file = open(energy_file_name, 'w')
    coord_file = open(coord_file_name, 'w')
    frc_file = open(frc_file_name, 'w')
    cell_file = open(cell_file_name, 'w')

    traj_coord_file_name = ''.join((sys_dir, '/active-pos-1.xyz'))
    traj_frc_file_name = ''.join((sys_dir, '/active-frc-1.xyz'))
    traj_cell_file_name = ''.join((sys_dir, '/active-1.cell'))
    traj_coord_file = open(traj_coord_file_name, 'w')
    traj_frc_file = open(traj_frc_file_name, 'w')
    traj_cell_file = open(traj_cell_file_name, 'w')

    for j in range(conv_iter):
      iter_dir = ''.join((work_dir, '/iter_', str(j)))
      data_dir = ''.join((iter_dir, '/03.cp2k_calc/sys_', str(i), '/data'))
      if ( j == 0 ):
        cmd = "cp type.raw %s" %(sys_dir)
        call.call_simple_shell(data_dir, cmd)
      energy_array, coord_array, frc_array, cell_array, virial_array = load_data.read_raw_data(data_dir)
      frames_num = len(energy_array)
      atoms_num = int(len(coord_array[0])/3)
      for k in range(frames_num):
        energy_file.write('%f\n' %(energy_array[k]))
        energy_cp2k.append(energy_array[k])
        frame_str = ''
        for l in range(len(coord_array[k])):
          if ( l == 0 ):
            frame_str = ''.join((frame_str, str(coord_array[k][l])))
          else:
            frame_str = ' '.join((frame_str, str(coord_array[k][l])))
        coord_file.write('%s\n' %(frame_str))

        frc_cp2k_k = []
        frc_x_cp2k_k = []
        frc_y_cp2k_k = []
        frc_z_cp2k_k = []

        frame_str = ''
        for l in range(len(frc_array[k])):
          if ( l == 0 ):
            frame_str = ''.join((frame_str, str(frc_array[k][l])))
          else:
            frame_str = ' '.join((frame_str, str(frc_array[k][l])))
          frc_cp2k_k.append(frc_array[k][l])
          if ( l%3 == 0 ):
            frc_x_cp2k_k.append(frc_array[k][l])
          elif ( l%3 == 1 ):
            frc_y_cp2k_k.append(frc_array[k][l])
          elif ( l%3 == 2 ):
            frc_z_cp2k_k.append(frc_array[k][l])
        frc_cp2k.append(frc_cp2k_k)
        frc_x_cp2k.append(frc_x_cp2k_k)
        frc_y_cp2k.append(frc_y_cp2k_k)
        frc_z_cp2k.append(frc_z_cp2k_k)
        frc_file.write('%s\n' %(frame_str))

        frame_str = ''
        for l in range(len(cell_array[k])):
          if ( l == 0 ):
            frame_str = ''.join((frame_str, str(cell_array[k][l])))
          else:
            frame_str = ' '.join((frame_str, str(cell_array[k][l])))
        cell_file.write('%s\n' %(frame_str))
    energy_file.close()
    coord_file.close()
    frc_file.close()
    cell_file.close()

    energy_array, coord_array, frc_array, cell_array, virial_array = load_data.read_raw_data(sys_dir)
    load_data.raw_data_to_set(1, False, sys_dir, energy_array, coord_array, frc_array, cell_array, virial_array)

    atoms = []
    type_raw = open(''.join((sys_dir, '/type.raw')), 'rb').read().split()
    for j in range(len(type_raw)):
      atoms.append(data_op.get_dic_keys(tot_atoms_type_dic, int(type_raw[j].decode())))

    traj_cell_file.write('#   Step   Time [fs]       Ax [Angstrom]       Ay [Angstrom]       Az [Angstrom]       Bx [Angstrom]       By [Angstrom]       Bz [Angstrom]       Cx [Angstrom]       Cy [Angstrom]       Cz [Angstrom]      Volume [Angstrom^3]\n')
    frames_num_tot = len(energy_array)
    for j in range(frames_num_tot):
      frc_array_j = frc_array[j].reshape(atoms_num, 3)
      coord_array_j = coord_array[j].reshape(atoms_num, 3)
      cell_array_j = cell_array[j].reshape(3,3)
      vol = np.linalg.det(cell_array_j)

      traj_coord_file.write('%8d\n' %(atoms_num))
      traj_coord_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', j, ', time =', j*0.5, ', E =', energy_array[j]))
      for k in range(atoms_num):
        traj_coord_file.write('%3s%21.10f%20.10f%20.10f\n' %(atoms[k], coord_array_j[k][0], coord_array_j[k][1], coord_array_j[k][2]))

      traj_frc_file.write('%8d\n' %(atoms_num))
      traj_frc_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', j, ', time =', j*0.5, ', E =', energy_array[j]))
      for k in range(atoms_num):
        traj_frc_file.write('%3s%21.10f%20.10f%20.10f\n' %(atoms[k], frc_array_j[k][0], frc_array_j[k][1], frc_array_j[k][2]))

      traj_cell_file.write('%8d%12.3f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%25.10f\n' \
                           %(j, j*0.5, cell_array_j[0][0], cell_array_j[0][1], cell_array_j[0][2], \
                             cell_array_j[1][0], cell_array_j[1][1], cell_array_j[1][2], \
                             cell_array_j[2][0], cell_array_j[2][1], cell_array_j[2][2], vol))

    traj_coord_file.close()
    traj_frc_file.close()
    traj_cell_file.close()

  str_print = 'Active data is written in %s' %(active_data_dir)
  print (data_op.str_wrap(str_print, 80), flush=True)

def run_iter(inp_file, deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, init_train_data, \
             init_data_num, work_dir, tot_atoms_type_dic, proc_num, host, device, usage):

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
    proc_num: int
      proc_num is the number of processors.
    host: 1-d string list
      host is the name of host of computational nodes.
    device: 2-d string list
      device is the gpu device name for each computational node.
    usage: 2-d float list
      usage is the memory usage of gpu devices for each computational node.
  Returns:
    none
  '''

  max_iter = model_devi_dic['max_iter']
  restart_iter = model_devi_dic['restart_iter']

  if ( restart_iter == 0 ):
    restart_data_num = init_data_num
  else:
    restart_data_num = model_devi_dic['restart_data_num']
  restart_stage = model_devi_dic['restart_stage']

  numb_test = deepmd_dic['training']['numb_test']
  model_type = deepmd_dic['training']['model_type']
  neuron = deepmd_dic['training']['neuron']
  shuffle_data = deepmd_dic['training']['shuffle_data']
  train_stress = deepmd_dic['training']['train_stress']

  model_devi_freq = int(model_devi_dic['model_devi_freq'])
  nsteps = int(lammps_dic['nsteps'])
  conv_new_data_num = int(nsteps/model_devi_freq*0.02)
  choose_new_data_num_limit = model_devi_dic['choose_new_data_num_limit']
  force_conv = model_devi_dic['force_conv']

  cp2k_exe = environ_dic['cp2k_exe']
  cp2k_env_file = environ_dic['cp2k_env_file']
  parallel_exe = environ_dic['parallel_exe']
  cuda_dir = environ_dic['cuda_dir']
  cp2k_job_num = environ_dic['cp2k_job_num']
  lmp_mpi_num = environ_dic['lmp_mpi_num']
  lmp_openmp_num = environ_dic['lmp_openmp_num']

  dp_path = sysinfo.get_dp_path(work_dir)
  lmp_path = sysinfo.get_lmp_path(work_dir)
  mpi_path = sysinfo.get_mpi_path(work_dir)

  data_num = []
  data_num.append(restart_data_num)

  np.random.seed(1234567890)

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
                                 descr_seed, fit_seed, tra_seed, neuron, model_type, sum(data_num))
      deepmd_run.run_deepmd(work_dir, i, parallel_exe, dp_path, host, device, usage, cuda_dir)
      check_deepff.write_restart_inp(inp_file, i, 1, sum(data_num), work_dir)

    if ( restart_stage == 0 or restart_stage == 1 ):
      #Perform lammps calculations
      print ('Step 2: lammps tasks', flush=True)

      lammps_run.gen_lmpmd_task(lammps_dic, work_dir, i, tot_atoms_type_dic)
      lammps_run.run_lmpmd(work_dir, i, lmp_path, mpi_path, lmp_mpi_num, lmp_openmp_num, device[0])
      check_deepff.write_restart_inp(inp_file, i, 2, sum(data_num), work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 ):
      #Perform lammps force calculations
      if ( restart_stage == 2 ):
        print ('Step 2: lammps tasks', flush=True)
      sys_num, atoms_type_multi_sys, atoms_num_tot = lammps_run.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      lammps_run.gen_lmpfrc_file(work_dir, i, atoms_num_tot, atoms_type_multi_sys)
      lammps_run.run_lmpfrc(work_dir, i, lmp_path, mpi_path, parallel_exe, proc_num, atoms_num_tot)
      check_deepff.write_restart_inp(inp_file, i, 3, sum(data_num), work_dir)

    if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 or restart_stage == 3 ):
      #Get force-force correlation and then choose new structures
      sys_num, atoms_type_multi_sys, atoms_num_tot = lammps_run.get_md_sys_info(lammps_dic, tot_atoms_type_dic)
      struct_index, success_ratio_sys, success_ratio, success_devi_ratio = \
      model_devi.choose_lmp_str(work_dir, i, atoms_type_multi_sys, force_conv)

      for j in range(len(success_ratio_sys)):
        print ('  The accurate ratio for system %d in iteration %d is %.2f%%' %(j, i, success_ratio_sys[j]*100), flush=True)

      print ('  The accurate ratio for whole %d systems in iteration %d is %.2f%% %.2f%%' \
             %(sys_num, i, success_ratio*100, success_devi_ratio*100), flush=True)

      if ( min(success_ratio_sys) >= 0.98 and success_ratio+success_devi_ratio >= 0.999 ):
        print (''.center(80,'*'), flush=True)
        print ('Cheers! deepff is converged!', flush=True)
        if ( i != 0 ):
          write_active_data(work_dir, i, tot_atoms_type_dic)
        exit()

      print ('Step 3: cp2k tasks', flush=True)
      #Perform cp2k calculation
      cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, i, atoms_type_multi_sys, atoms_num_tot, \
                             struct_index, conv_new_data_num, choose_new_data_num_limit, train_stress)
      cp2k_run.run_cp2kfrc(work_dir, i, cp2k_exe, parallel_exe, cp2k_env_file, cp2k_job_num, proc_num, atoms_num_tot)

      #Dump new data of cp2k
      for j in range(sys_num):
        file_dir = ''.join((work_dir, '/iter_', str(i), '/03.cp2k_calc/sys_', str(j)))
        load_data.load_data_from_sepfile(file_dir, 'task_', 'cp2k', tot_atoms_type_dic)
        cp2k_data_dir = ''.join((file_dir, '/data'))
        if ( os.path.exists(cp2k_data_dir) ):
          energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(cp2k_data_dir)
          train_data_num = load_data.raw_data_to_set(1, shuffle_data, cp2k_data_dir, energy_array, coord_array, frc_array, box_array, virial_array)
          data_num.append(train_data_num)

      print ('  Success: dump new raw data of cp2k', flush=True)
      check_deepff.write_restart_inp(inp_file, i+1, 0, sum(data_num), work_dir)
      restart_stage = 0

    if ( i == max_iter-1 ):
      log_info.log_error('Active learning does not converge')
      write_active_data(work_dir, i+1, tot_atoms_type_dic)

def kernel(work_dir, inp_file):

  '''
  kernel: kernel function to do active learning.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    inp_file: string
      inp_file is the deepff input file
  '''

  proc_num, host, ssh = sysinfo.get_host(work_dir)
  device, usage = sysinfo.analyze_gpu(host, ssh, work_dir)

  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']

  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  dump_input(work_dir, inp_file, deepff_key)

  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  restart_iter = model_devi_dic['restart_iter']
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

  init_train_data, init_data_num = dump_init_data(work_dir, deepmd_dic, restart_iter, train_stress, tot_atoms_type_dic)

  print ('Check input file: no error in %s' %(inp_file), flush=True)
  print ('Initial training data:', flush=True)
  for i in range(len(init_train_data)):
    print ('%s' %(data_op.str_wrap(init_train_data[i], 80)), flush=True)

  run_iter(inp_file, deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, init_train_data, \
           init_data_num, work_dir, tot_atoms_type_dic, proc_num, host, device, usage)

if __name__ == '__main__':

  from CP2K_kit.deepff import active
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  inp_file = 'input.inp'
  max_cycle = 100
  active.kernel(work_dir, inp_file)
