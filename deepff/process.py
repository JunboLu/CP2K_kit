#! /usr/env/bin python

import os
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import file_tools
from CP2K_kit.tools import read_input
from CP2K_kit.tools import traj_info
from CP2K_kit.deepff import load_data
from CP2K_kit.deepff import gen_lammps_task

def get_sys_num(exe_dir):

  '''
  get_sys_num: get the number of systems

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    sys_num: int
      sys_num is the number of systems.
  '''

  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(exe_dir, cmd))

  return sys_num

def get_data_num(exe_dir):

  '''
  get_data_num: get the number of data

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    data_num: int
      data_num is the number of data.
  '''

  cmd = "ls | grep %s" % ('data_')
  data_num = len(call.call_returns_shell(exe_dir, cmd))

  return data_num

def get_task_num(exe_dir, get_task_dir=False):

  '''
  get_task_num: get the number of tasks in a system

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    task_num: int
      task_num is the number of tasks in a system.
  '''

  cmd = "ls | grep %s" % ('task_')
  task_dir = call.call_returns_shell(exe_dir, cmd)
  task_num = len(task_dir)

  if get_task_dir:
    return task_num, task_dir
  else:
    return task_num

def get_lmp_model_num(exe_dir):

  '''
  get_lmp_model_num: get the number of models in lammps directory.

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    model_num: int
      model_num is the number of models in lammps directory.
  '''

  cmd = "ls | grep %s" % ("'model_[0-9]'")
  model_num = len(call.call_returns_shell(exe_dir, cmd))

  return model_num

def get_deepmd_model_num(exe_dir):

  '''
  get_deepmd_model_num: get the number of models in deepmd directory.

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    model_num: int
      model_num is the number of models in deepmd directory.
  '''

  model_num = len(call.call_returns_shell(exe_dir, "ls -ll |awk '/^d/ {print $NF}'"))

  return model_num

def get_traj_num(exe_dir):

  '''
  get_traj_num: get the number of frames

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
  Returns:
    traj_num: int
      traj_num is the number of frames.
  '''

  cmd = "ls | grep %s" % ('traj_')
  traj_num = len(call.call_returns_shell(exe_dir, cmd))

  return traj_num

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
  active_learn_dic = job_type_param[3]
  environ_dic = job_type_param[4]

  return deepmd_dic, lammps_dic, cp2k_dic, active_learn_dic, environ_dic

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
      traj_type = train_dic[key]['traj_type']
      if ( traj_type == 'md' ):
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
      elif ( traj_type == 'mtd' ):
        data_dir = train_dic[key]['data_dir']
        task_dir_prefix = train_dic[key]['task_dir_prefix']
        start_frame = train_dic[key]['start_frame']
        coord_file = ''.join((data_dir, '/', task_dir_prefix, str(start_frame), '/coord'))
        atoms_num = len(open(coord_file).readlines())
        for i in range(atoms_num):
          line_i = linecache.getline(coord_file, i+1)
          line_i_split = data_op.split_str(line_i, ' ', '\n')
          atoms.append(line_i_split[0])
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
          choosed_index = data_op.gen_list(start, end, 1)
          choosed_index_array = np.array(choosed_index)
          np.random.shuffle(choosed_index_array)
          choosed_index = list(choosed_index_array[0:choosed_num])
          load_data.load_data_from_sepfile(data_dir, save_dir, task_dir_prefix, proj_name, tot_atoms_type_dic, \
                                           sorted(choosed_index), out_file_name)
        energy_array, coord_array, frc_array, box_array, virial_array = load_data.read_raw_data(save_dir)
        init_data_num_part, init_test_data_num_part = load_data.raw_data_to_set(parts, shuffle_data, save_dir, energy_array, \
                                                                                coord_array, frc_array, box_array, virial_array)
      init_data_num = init_data_num+init_data_num_part
      i = i+1

  if ( 'set_data_dir' in train_dic.keys() ):
    init_train_data.append(os.path.abspath(os.path.expanduser(train_dic['set_data_dir'])))
    energy_npy_file = ''.join((os.path.abspath(os.path.expanduser(train_dic['set_data_dir'])), '/set.000/energy.npy'))
    set_data_num = len(np.load(energy_npy_file))
    init_data_num = init_data_num+set_data_num

  return init_train_data, init_data_num

def check_deepff_run(work_dir, iter_id):

  '''
  check_deepff_run: check the running state of deepff

  Args:
   work_dir: string
      work_dir is workding directory.
    iter_id: int
      iter_id is current iteration number.
  Returns:
    failure_model: 1-d int list
      failure_model is the id of failure models.
  '''

  train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  model_num = get_deepmd_model_num(train_dir)
  failure_model = []
  for i in range(model_num):
    model_dir = ''.join((train_dir, '/', str(i)))
    lcurve_file = ''.join((model_dir, '/lcurve.out'))
    whole_line_num = len(open(lcurve_file).readlines())
    choosed_line_num = int(0.1*whole_line_num)
    start_line = whole_line_num-choosed_line_num
    force_trn = []
    for j in range(choosed_line_num):
      line = linecache.getline(lcurve_file, start_line+j+1)
      line_split = data_op.split_str(line, ' ')
      if ( data_op.eval_str(line_split[0]) == 1 and len(line_split) >= 8 ):
        force_trn.append(float(line_split[6]))
    linecache.clearcache()
    force_max = max(force_trn)
    force_min = min(force_trn)
    force_avg = np.mean(np.array(force_trn))
    if ( ((force_max-force_min) >= 0.05 and force_max >= 0.09) or force_avg >= 0.09 ):
      failure_model.append(i)

  return failure_model

def get_md_sys_info(lmp_dic, tot_atoms_type_dic):

  '''
  get_md_sys_info: get the system information for lammps md.

  Args:
    lmp_dic: dictionary
      lmp_dic contains parameters for lammps.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    sys_num: int
      sys_num is the number of systems.
    atoms_type_multi_sys: 2-d dictionary, dim = (num of lammps systems) * (num of atom types)
      atoms_type_multi_sys is the atoms type for multi-systems.
      example: {0:{'O':1,'H':2,'N':3},1:{'O':1,'S':2,'N':3}}
    atoms_num_tot: dictionary
      atoms_num_tot contains number of atoms for different systems.
      Example: {0:3, 1:3}
    use_mtd_tot: bool
      use_mtd_tot is whethet using metadynamics for whole systems.
  '''

  atoms_type_multi_sys = []
  atoms_num_tot = []
  use_mtd_tot = []
  sys_num = 0
  for key in lmp_dic:
    if 'system' in key:
      sys_num = sys_num + 1

  for i in range(sys_num):
    sys = 'system' + str(i)
    box_file_name = lmp_dic[sys]['box']
    coord_file_name = lmp_dic[sys]['coord']
    use_mtd = lmp_dic[sys]['use_mtd']
    tri_cell_vec, atoms, x, y, z = gen_lammps_task.get_box_coord(box_file_name, coord_file_name)
    atoms_type = data_op.list_replicate(atoms)
    atoms_type_dic = OrderedDict()
    for j in atoms_type:
      if j in tot_atoms_type_dic.keys():
        atoms_type_dic[j] = tot_atoms_type_dic[j]+1
      else:
        log_info.log_error('Input error: %s atom type in system %d is not trained, please check deepff/lammps/system' %(j, i))
        exit()
    atoms_type_multi_sys.append(atoms_type_dic)
    atoms_num_tot.append(len(atoms))
    use_mtd_tot.append(use_mtd)

  return sys_num, atoms_type_multi_sys, atoms_num_tot, use_mtd_tot

if __name__ == '__main__':
  from CP2K_kit.deepff import process
  check_deepff_run('/home/lujunbo/WORK/Deepmd/C2H6/TRAINING/train_md_mtd/active_mtd', 0)
