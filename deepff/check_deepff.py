#! /usr/env/bin python

import os
import linecache
import multiprocessing
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import file_tools

def check_inp(deepmd_dic, lmp_dic, cp2k_dic, model_devi_dic, environ_dic):

  '''
  check_inp: check the deepff input file

  Args:
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
  Returns:
    deepmd_dic: dictionary
      deepmd_dic is the revised deepmd_dic.
    lmp_dic: dictionary
      lmp_dic is the revised lmp_dic.
    cp2k_dic: dictionary
      cp2k_dic is the revised cp2k_dic.
    model_devi_dic: dictionary
      model_devi_dic is the revised model_devi_dic.
    environ_dic: dictionary
      environ_dic is the revised environ_dic.
  '''

  #Check parameters for deepmd-kit
  if ( 'model' not in deepmd_dic.keys() ):
    log_info.log_error('Input error: no model, please set deepff/deepmd/model')
    exit()
  else:
    if ( 'type_map' in deepmd_dic['model'].keys() ):
      if ( all(data_op.eval_str(i) == 0 for i in deepmd_dic['model']['type_map']) ):
        pass
      else:
        log_info.log_error('Input error: type_map should be string, please check or reset deepff/deepmd/model/type_map')
        exit()
    else:
      log_info.log_error('Input error: no type_map, please set deepff/deepmd/model/type_map')
      exit()
    if ( 'descriptor' not in deepmd_dic['model'].keys() ):
      log_info.log_error('Input error: no descriptor, please set deepff/deepmd/model/descriptor')
      exit()
    else:
      valid_type = ['se_a', 'se_r', 'se_ar']
      if ( 'type' in deepmd_dic['model']['descriptor'].keys() ):
        descr_type = deepmd_dic['model']['descriptor']['type']
        if ( descr_type in valid_type ):
          pass
        else:
          log_info.log_error('Input error: %s is not supported for deepff/deepmd/model/descriptor/type')
          exit()

      if ( 'sel' in deepmd_dic['model']['descriptor'].keys() ):
        sel = deepmd_dic['model']['descriptor']['sel']
        if ( all(data_op.eval_str(i) == 1 for i in sel) ):
          deepmd_dic['model']['descriptor']['sel'] = [int(x) for x in sel]
        else:
          log_info.log_error('Input error: sel shoule be list of integer, please check or reset deepff/deepmd/model/descriptor/sel')
          exit()
      else:
        log_info.log_error('Input error: no sel, please set deepff/deepmd/model/descriptor/sel')
        exit()

      if ( 'rcut_smth' in deepmd_dic['model']['descriptor'].keys() ):
        rcut_smth = deepmd_dic['model']['descriptor']['rcut_smth']
        if ( data_op.eval_str(rcut_smth) == 1 or data_op.eval_str(rcut_smth) == 2 ):
          deepmd_dic['model']['descriptor']['rcut_smth'] = float(rcut_smth)
        else:
          log_info.log_error('Input error: rcut_smth shoule be float, please check or reset deepff/deepmd/model/descriptor/rcut_smth')
          exit()
      else:
        deepmd_dic['model']['descriptor']['rcut_smth'] = 5.0

      if ( 'rcut' in deepmd_dic['model']['descriptor'].keys() ):
        rcut = deepmd_dic['model']['descriptor']['rcut']
        if ( data_op.eval_str(rcut) == 1 or data_op.eval_str(rcut) == 2 ):
          deepmd_dic['model']['descriptor']['rcut'] = float(rcut)
        else:
          log_info.log_error('Input error: rcut should be float, please check or reset deepff/deepmd/model/descriptor/rcut')
          exit()
      else:
        deepmd_dic['model']['descriptor']['rcut'] = 6.0

      if ( 'neuron' in deepmd_dic['model']['descriptor'].keys() ):
        neuron_encode = deepmd_dic['model']['descriptor']['neuron']
        if ( all(data_op.eval_str(i) == 1 for i in neuron_encode) ):
          deepmd_dic['model']['descriptor']['neuron'] = [int(x) for x in neuron_encode]
        else:
          log_info.log_error('Input error: neuron error, please check deepff/deepmd/model/descriptor/neuron')
          exit()
      else:
        deepmd_dic['model']['descriptor']['neuron'] = [25, 50, 100]

      deepmd_dic['model']['descriptor']['resnet_dt'] = False

      if ( 'axis_neuron' in deepmd_dic['model']['descriptor'].keys() ):
        axis_neuron = deepmd_dic['model']['descriptor']['axis_neuron']
        if ( data_op.eval_str(axis_neuron) == 1 ):
          deepmd_dic['model']['descriptor']['axis_neuron'] = int(axis_neuron)
        else:
          log_info.log_error('Input error: axis_neuron should be list of integer, please check deepff/deepmd/model/descriptor/axis_neuron')
          exit()
      else:
        deepmd_dic['model']['descriptor']['axis_neuron'] = 16

    deepmd_dic['model']['fitting_net'] = OrderedDict()
    deepmd_dic['model']['fitting_net']['resnet_dt'] = True

  if ( 'learning_rate' not in deepmd_dic.keys() ):
    log_info.log_error('Input error: no learning_rate, please set deepff/deepmd/learning_rate')
    exit()
  else:
    if ( 'type' in deepmd_dic['learning_rate'].keys() ):
      decay_type = deepmd_dic['learning_rate']['type']
      if ( data_op.eval_str(decay_type) == 0 ):
        pass
      else:
        log_info.log_error('Input error: type in learning_rate should be string, please check or reset deepff/deepmd/learning_rate/type')
        exit()
    else:
      deepmd_dic['learning_rate']['type'] = 'exp'

    if ( 'start_lr' in deepmd_dic['learning_rate'].keys() ):
      start_lr = deepmd_dic['learning_rate']['start_lr']
      if ( data_op.eval_str(start_lr) == 2 ):
        deepmd_dic['learning_rate']['start_lr'] = float(start_lr)
      else:
        log_info.log_error('Input error: start_lr should be float, please check deepff/deepmd/learning_rate/start_lr')
        exit()
    else:
      deepmd_dic['learning_rate']['start_lr'] = 0.001

    if ( 'stop_lr' in deepmd_dic['learning_rate'].keys() ):
      stop_lr = deepmd_dic['learning_rate']['stop_lr']
      if ( data_op.eval_str(stop_lr) == 2 ):
        deepmd_dic['learning_rate']['stop_lr'] = float(stop_lr)
      else:
        log_info.log_error('Input error: stop_lr should be float, please check deepff/deepmd/learning_rate/stop_lr')
        exit()
    else:
      deepmd_dic['learning_rate']['stop_lr'] = 1e-8

  if ( 'loss' not in deepmd_dic.keys() ):
    log_info.log_error('Input error: no loss, please check or set deepff/deepmd/loss')
    exit()
  else:
    if ( 'start_pref_e' in deepmd_dic['loss'].keys() ):
      start_pref_e = deepmd_dic['loss']['start_pref_e']
      if ( data_op.eval_str(start_pref_e) == 1 or data_op.eval_str(start_pref_e) == 2 ):
        deepmd_dic['loss']['start_pref_e'] = float(start_pref_e)
      else:
        log_info.log_error('Input error: start_pref_e should be float, please check deepff/deepmd/loss/start_pref_e')
        exit()
    else:
      deepmd_dic['loss']['start_pref_e'] = 0.02

    if ( 'limit_pref_e' in deepmd_dic['loss'].keys() ):
      limit_pref_e = deepmd_dic['loss']['limit_pref_e']
      if ( data_op.eval_str(limit_pref_e) == 1 or data_op.eval_str(limit_pref_e) == 2 ):
        deepmd_dic['loss']['limit_pref_e'] = float(limit_pref_e)
      else:
        log_info.log_error('Input error: limit_pref_e should be float, please check deepff/deepmd/loss/limit_pref_e')
        exit()
    else:
      deepmd_dic['loss']['limit_pref_e'] = 1.0

    if ( 'start_pref_f' in deepmd_dic['loss'].keys() ):
      start_pref_f = deepmd_dic['loss']['start_pref_f']
      if ( data_op.eval_str(start_pref_f) == 1 or data_op.eval_str(start_pref_f) == 2 ):
        deepmd_dic['loss']['start_pref_f'] = float(start_pref_f)
      else:
        log_info.log_error('Input error: start_pref_f should be float, please check deepff/deepmd/loss/start_pref_f')
        exit()
    else:
      deepmd_dic['loss']['start_pref_f'] = 1000.0

    if ( 'limit_pref_f' in deepmd_dic['loss'].keys() ):
      limit_pref_f = deepmd_dic['loss']['limit_pref_f']
      if ( data_op.eval_str(limit_pref_f) == 1 or data_op.eval_str(limit_pref_f) == 2 ):
        deepmd_dic['loss']['limit_pref_f'] = float(limit_pref_f)
      else:
        log_info.log_error('Input error: limit_pref_f should be float, please check deepff/deepmd/loss/limit_pref_f')
        exit()
    else:
      deepmd_dic['loss']['limit_pref_f'] = 1.0

    if ( 'start_pref_v' in deepmd_dic['loss'].keys() ):
      start_pref_v = deepmd_dic['loss']['start_pref_v']
      if ( data_op.eval_str(start_pref_v) == 1 or data_op.eval_str(start_pref_v) == 2 ):
        deepmd_dic['loss']['start_pref_v'] = float(start_pref_v)
      else:
        log_info.log_error('Input error: start_pref_v should be float, please check deepff/deepmd/loss/start_pref_v')
        exit()
    else:
      deepmd_dic['loss']['start_pref_v'] = 0.0

    if ( 'limit_pref_v' in deepmd_dic['loss'].keys() ):
      limit_pref_v = deepmd_dic['loss']['limit_pref_v']
      if ( data_op.eval_str(limit_pref_v) == 1 or data_op.eval_str(limit_pref_v) == 2 ):
        deepmd_dic['loss']['limit_pref_v'] = float(limit_pref_v)
      else:
        log_info.log_error('Input error: limit_pref_v should, please check deepff/deepmd/loss/limit_pref_v')
        exit()
    else:
      deepmd_dic['loss']['limit_pref_v'] = 0.0

  if ( 'training' not in deepmd_dic.keys() ):
    log_info.log_error('Input error: no training found, please check or set deepff/deepmd/training')
    exit()
  else:
    for i in deepmd_dic['training'].keys():
      if ( 'system' in i ):
        if ( 'traj_coord_file' in deepmd_dic['training'][i] ):
          traj_coord_file = deepmd_dic['training'][i]['traj_coord_file']
          if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
            deepmd_dic['training'][i]['traj_coord_file'] = os.path.abspath(traj_coord_file)
          else:
            log_info.log_error('Input error: %s does not exist, please check deepff/deepmd/training/system/traj_coord_file' %(traj_coord_file))
            exit()
        else:
          log_info.log_error('Input error: no coordination trajectory file, please check deepff/deepmd/training/system/traj_coord_file')
          exit()
        if ( 'traj_frc_file' in deepmd_dic['training'][i] ):
          traj_frc_file = deepmd_dic['training'][i]['traj_frc_file']
          if ( os.path.exists(os.path.abspath(traj_frc_file)) ):
            deepmd_dic['training'][i]['traj_frc_file'] = os.path.abspath(traj_frc_file)
          else:
            log_info.log_error('Input error: %s does not exist, please check deepff/deepmd/training/system/traj_frc_file' %(traj_frc_file))
            exit()
        else:
          log_info.log_error('Input error: no force trajectory file, please check deepff/deepmd/training/system/traj_frc_file')
          exit()
        line_num = file_tools.grep_line_num("'PDB file'", traj_coord_file, os.getcwd())
        if ( line_num == 0 ):
          coord_file_type = 'coord_xyz'
        else:
          coord_file_type = 'coord_pdb'
        if ( coord_file_type == 'coord_xyz' ):
          if ( 'traj_cell_file' in deepmd_dic['training'][i] ):
            traj_cell_file = deepmd_dic['training'][i]['traj_cell_file']
            if ( os.path.exists(os.path.abspath(traj_cell_file)) ):
              deepmd_dic['training'][i]['traj_cell_file'] = os.path.abspath(traj_cell_file)
            else:
              log_info.log_error('Input error: %s does not exist, please check deepff/deepmd/training/system/traj_cell_file' %(traj_cell_file))
              exit()
          else:
            log_info.log_error('Input error: no cell trajectory file, please check deepff/deepmd/training/system/traj_cell_file')
            exit()
        else:
          deepmd_dic['training'][i]['traj_cell_file'] = 'none'
        if ( 'traj_stress_file' in deepmd_dic['training'][i] ):
          traj_stress_file = deepmd_dic['training'][i]['traj_stress_file']
          if ( os.path.exists(os.path.abspath(traj_stress_file)) ):
            deepmd_dic['training'][i]['traj_stress_file'] = os.path.abspath(traj_stress_file)
          else:
            log_info.log_error('Input error: %s does not exist, please check deepff/deepmd/training/system/traj_stress_file' %(traj_stress_file))
            exit()
        else:
          deepmd_dic['training'][i]['traj_stress_file'] = 'none'
        atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_id, end_id, time_step = \
        traj_info.get_traj_info(traj_coord_file, coord_file_type)
        if ( 'start_frame' in deepmd_dic['training'][i] ):
          start_frame = deepmd_dic['training'][i]['start_frame']
          if ( data_op.eval_str(start_frame) == 1 ):
            deepmd_dic['training'][i]['start_frame'] = int(start_frame)
          else:
            log_info.log_error('Input error: start_frame should be integer, please check deepff/deepmd/training/system/start_frame')
            exit()
        else:
          deepmd_dic['training'][i]['start_frame'] = start_id
        if ( 'end_frame' in deepmd_dic['training'][i] ):
          end_frame = deepmd_dic['training'][i]['end_frame']
          if ( data_op.eval_str(end_frame) == 1 ):
            deepmd_dic['training'][i]['end_frame'] = int(end_frame)
          else:
            log_info.log_error('Input error: end_frame should be integer, please check deepff/deepmd/training/system/end_frame')
            exit()
        else:
          deepmd_dic['training'][i]['end_frame'] = end_id
        if ( 'choosed_frame_num' in deepmd_dic['training'][i] ):
          choosed_frame_num = deepmd_dic['training'][i]['choosed_frame_num']
          if ( data_op.eval_str(choosed_frame_num) == 1 ):
            deepmd_dic['training'][i]['choosed_frame_num'] = int(choosed_frame_num)
          else:
            log_info.log_error('Input error: choosed_frame_num should be integer, please check deepff/deepmd/training/system/choosed_frame_num')
            exit()
        else:
          deepmd_dic['training'][i]['choosed_frame_num'] = int((end_id-start_id)/each)+1
        if ( 'set_parts' in deepmd_dic['training'][i] ):
          set_parts = deepmd_dic['training'][i]['set_parts']
          if ( data_op.eval_str(set_parts) == 1 ):
            deepmd_dic['training'][i]['set_parts'] = int(set_parts)
          else:
            log_info.log_error('Input error: set_parts should be integer, please check deepff/deepmd/training/system/set_parts')
            exit()
        else:
          deepmd_dic['training'][i]['set_parts'] = 1

    if ( 'use_prev_model' in deepmd_dic['training'].keys() ):
      use_prev_model = deepmd_dic['training']['use_prev_model']
      use_prev_model_bool = data_op.str_to_bool(use_prev_model)
      if ( isinstance(use_prev_model_bool, bool) ):
        deepmd_dic['training']['use_prev_model'] = use_prev_model_bool
      else:
        log_info.log_error('Input error: use_prev_model should be bool, please check or reset deepff/deepmd/training/use_prev_model')
        exit()
    else:
      deepmd_dic['training']['use_prev_model'] = False

    if ( 'shuffle_data' in deepmd_dic['training'].keys() ):
      shuffle_data = deepmd_dic['training']['shuffle_data']
      shuffle_data_bool = data_op.str_to_bool(shuffle_data)
      if ( isinstance(shuffle_data_bool, bool) ):
        deepmd_dic['training']['shuffle_data'] = shuffle_data_bool
      else:
        log_info.log_error('Input error: shuffle_data should be bool, please check or reset deepff/deepmd/training/shuffle_data')
        exit()
    else:
      deepmd_dic['training']['shuffle_data'] = False

    if ( 'train_stress' in deepmd_dic['training'].keys() ):
      train_stress = deepmd_dic['training']['train_stress']
      train_stress_bool = data_op.str_to_bool(train_stress)
      if ( isinstance(train_stress_bool, bool) ):
        deepmd_dic['training']['train_stress'] = train_stress_bool
      else:
        log_info.log_error('Input error: train_stress should be bool, please check or reset deepff/deepmd/training/train_stress')
        exit()
    else:
      deepmd_dic['training']['train_stress'] = False

    if ( 'set_data_dir' in deepmd_dic['training'].keys() ):
      set_data_dir = deepmd_dic['training']['set_data_dir']
      if ( os.path.exists(os.path.abspath(set_data_dir)) ):
        deepmd_dic['training']['set_data_dir'] = os.path.abspath(set_data_dir)
      else:
        log_info.log_error('Input error: %s directory does not exist, please check or reset deepff/deepmd/training/set_data_dir')
        exit()

    if ( 'model_type' in deepmd_dic['training'].keys() ):
      model_type = deepmd_dic['training']['model_type']
      if ( model_type == 'use_seed' or model_type == 'use_node' ):
        pass
      else:
        log_info.log_error('Input error: only use_seed and use_node are supported for model_type, please check deepff/deepmd/training/model_type')
        exit()
    else:
      deepmd_dic['training']['model_type'] = 'use_seed'

    if ( deepmd_dic['training']['model_type'] == 'use_seed' ):
      if ( 'seed_num' in deepmd_dic['training'].keys() ):
        seed_num = deepmd_dic['training']['seed_num']
        if ( data_op.eval_str(seed_num) == 1 ):
          deepmd_dic['training']['seed_num'] = int(seed_num)
        else:
          log_info.log_error('Input error: seed_num should be integer, please check or reset deepff/deepmd/training/seed_num')
          exit()
      else:
        deepmd_dic['training']['seed_num'] = 2

    if ( 'neuron' in deepmd_dic['training'].keys() ):
      neuron_list = deepmd_dic['training']['neuron']
      if ( deepmd_dic['training']['model_type'] == 'use_node' ):
        neuron = []
        tmp_str = data_op.comb_list_2_str(neuron_list, ' ')
        tmp_list = data_op.split_str(tmp_str, '...')

        for i in range(len(tmp_list)):
          neuron_i = data_op.split_str(tmp_list[i], ' ')
          if ( all(data_op.eval_str(j) == 1 for j in neuron_i) ):
            neuron.append([int(x) for x in neuron_i])
          else:
            log_info.log_error('Input error: neuron should be list of integer, please check or reset deepff/deepmd/training/neuron')
            exit()
      elif ( deepmd_dic['training']['model_type'] == 'use_seed' ):
        if ( all(data_op.eval_str(j) == 1 for j in neuron_list) ):
          neuron = [int(x) for x in neuron_list]
        else:
          log_info.log_error('Input error: neuron should be list of integer, please check or reset deepff/deepmd/training/neuron')
          exit()
      deepmd_dic['training']['neuron'] = neuron
    else:
      log_info.log_error('Input error: no neuron, please set deepff/deepmd/training/neuron')
      exit()

    if ( 'fix_stop_batch' in deepmd_dic['training'].keys() ):
      fix_stop_batch = deepmd_dic['training']['fix_stop_batch']
      fix_stop_batch_bool = data_op.str_to_bool(fix_stop_batch)
      if ( isinstance(fix_stop_batch_bool, bool) ):
        deepmd_dic['training']['fix_stop_batch'] = fix_stop_batch_bool
      else:
        log_info.log_error('Input error: fix_stop_batch should be bool, please check or reset deepff/deepmd/training/fix_stop_batch')
        exit()
    else:
      deepmd_dic['training']['fix_stop_batch'] = False

    fix_stop_batch = deepmd_dic['training']['fix_stop_batch']
    if fix_stop_batch:
      if ( 'decay_steps' in deepmd_dic['learning_rate'].keys() ):
        decay_steps = deepmd_dic['learning_rate']['decay_steps']
        if ( data_op.eval_str(decay_steps) == 1 ):
          deepmd_dic['learning_rate']['decay_steps'] = int(decay_steps)
        else:
          log_info.log_error('Input error: the decay_steps should be integer, please check or reset deepff/deepmd/learning_rate/decay_steps')
          exit()
      else:
        deepmd_dic['learning_rate']['decay_steps'] = 5000
      if ( 'stop_batch' in deepmd_dic['training'].keys() ):
        stop_batch = deepmd_dic['training']['stop_batch']
        if ( data_op.eval_str(stop_batch) == 1 ):
          deepmd_dic['training']['stop_batch'] = int(stop_batch)
        else:
          log_info.log_error('Input error: the stop_batch should be integer, please check or reset deepff/deepmd/training/stop_batch')
          exit()
      else:
        deepmd_dic['training']['stop_batch'] = 1000000
    else:
      if ( 'epoch_num' in deepmd_dic['training'].keys() ):
        epoch_num = deepmd_dic['training']['epoch_num']
        if ( data_op.eval_str(epoch_num) == 1 ):
          deepmd_dic['training']['epoch_num'] = int(epoch_num)
        else:
          log_info.log_error('Input error: the number of epoch should be integer, please check or reset deepff/deepmd/training/epoch_num')
          exit()
      else:
        deepmd_dic['training']['epoch_num'] = 200

    if ( 'batch_size' in deepmd_dic['training'].keys() ):
      batch_size = deepmd_dic['training']['batch_size']
      if ( data_op.eval_str(batch_size) == 1 ):
        deepmd_dic['training']['batch_size'] = int(batch_size)
      else:
        log_info.log_error('Input error: batch_size shoule be integer, please check or reset deepff/deepmd/training/batch_size')
        exit()
    else:
      deepmd_dic['training']['batch_size'] = 1

    if ( 'disp_freq' in deepmd_dic['training'].keys() ):
      disp_freq = deepmd_dic['training']['disp_freq']
      if ( data_op.eval_str(disp_freq) == 1 ):
        deepmd_dic['training']['disp_freq'] = int(disp_freq)
      else:
        log_info.log_error('Input error: disp_freq should be integer, please check or reset deepff/deepmd/training/disp_freq')
        exit()
    else:
      deepmd_dic['training']['disp_freq'] = 100

    if ( 'numb_test' in deepmd_dic['training'].keys() ):
      numb_test = deepmd_dic['training']['numb_test']
      if ( data_op.eval_str(numb_test) == 1 ):
        deepmd_dic['training']['numb_test'] = int(numb_test)
      else:
        log_info.log_error('Input error: numb_test should be integer, please check or reset deepff/deepmd/training/numb_test')
        exit()
    else:
      deepmd_dic['training']['numb_test'] = 10

    if ( 'save_freq' in deepmd_dic['training'].keys() ):
      save_freq = deepmd_dic['training']['save_freq']
      if ( data_op.eval_str(save_freq) == 1 ):
        deepmd_dic['training']['save_freq'] = int(save_freq)
      else:
        log_info.log_error('Input error: save_freq should be integer, please check or reset deepff/deepmd/training/save_freq')
        exit()
    else:
      deepmd_dic['training']['save_freq'] = 1000

    if ( 'disp_training' in deepmd_dic['training'].keys() ):
      disp_training = deepmd_dic['training']['disp_training']
      disp_training_bool = data_op.str_to_bool(disp_training)
      if ( isinstance(disp_training_bool, bool) ):
        deepmd_dic['training']['disp_training'] = disp_training_bool
      else:
        log_info.log_error('Input error: disp_training should be bool, please check or reset deepff/deepmd/training/disp_training')
        exit()
    else:
      deepmd_dic['training']['disp_training'] = True

    deepmd_dic['training']['set_prefix'] = 'set'
    deepmd_dic['training']['disp_file'] = 'lcurve.out'
    deepmd_dic['training']['load_ckpt'] = 'model.ckpt'
    deepmd_dic['training']['save_ckpt'] = 'model.ckpt'
    deepmd_dic['training']['time_training'] = True
    deepmd_dic['training']['profiling'] = False
    deepmd_dic['training']['profiling_file'] = 'timeline.json'

  #Check parameters for lammps
  if ( 'nsteps' in lmp_dic.keys() ):
    nsteps = lmp_dic['nsteps']
    if ( data_op.eval_str(nsteps) == 1 ):
      if ( int(nsteps) < 1000 ):
        lmp_dic['nsteps'] = '1000'
    else:
      log_info.log_error('Input error: nsteps should be integer, please check or reset deepff/lammps/nsteps')
      exit()
  else:
    lmp_dic['nsteps'] = '10000'

  if ( 'write_restart_freq' in lmp_dic.keys() ):
    write_restart_freq = lmp_dic['write_restart_freq']
    if ( data_op.eval_str(write_restart_freq) == 1 ):
      pass
    else:
      log_info.log_error('Input error: write_restart_freq should be integer, please check or reset deepff/lammps/write_restart_freq')
      exit()
  else:
    lmp_dic['write_restart_freq'] = '1000'

  if ( 'model_devi_freq' in model_devi_dic.keys() ):
    model_devi_freq = model_devi_dic['model_devi_freq']
    if ( data_op.eval_str(model_devi_freq) == 1 ):
      pass
    else:
      log_info.log_error('Input error: model_devi_freq should be integer, please check or reset deepff/lammps/model_devi_freq')
      exit()
  else:
    model_devi_dic['model_devi_freq'] = '10'

  model_devi_freq = model_devi_dic['model_devi_freq']

  lmp_dic['thermo_freq'] = model_devi_freq
  lmp_dic['dump_freq'] = model_devi_freq

  if ( 'time_step' in lmp_dic.keys() ):
    time_step = lmp_dic['time_step']
    if ( data_op.eval_str(time_step) == 1 or data_op.eval_str(time_step) == 2 ):
      pass
    else:
      log_info.log_error('Input error: time_step should be integer or float, please check or reset deepff/lammps/time_step')
  else:
    lmp_dic['time_step'] = '0.0005'

  if ( 'tau_t' in lmp_dic.keys() ):
    tau_t = lmp_dic['tau_t']
    if ( data_op.eval_str(tau_t) == 1 or data_op.eval_str(tau_t) == 2 ):
      pass
    else:
      log_info.log_error('Input error: tau_t should be integer or float, please check or reset deepff/lammps/tau_t')
      exit()
  else:
    lmp_dic['tau_t'] = '%f' %(float(lmp_dic['time_step'])*200)

  if ( 'tau_p' in lmp_dic.keys() ):
    tau_p = lmp_dic['tau_p']
    if ( data_op.eval_str(tau_p) == 2 ):
      pass
    else:
      log_info.log_error('Input error: tau_p should be integer or float, please check or reset deepff/lammps/tau_p')
      exit()
  else:
    lmp_dic['tau_p'] = '%f' %(float(lmp_dic['time_step'])*200)

  if ( 'temp' in lmp_dic.keys() ):
    temp = lmp_dic['temp']
    if ( isinstance(temp, list) ):
      if ( all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in temp)):
        lmp_dic['temp'] = [float(i) for i in temp]
      else:
        log_info.log_error('Input error: multipole temperature should be list of float, please check or reset deepff/lammps/temp')
        exit()
    else:
      if ( data_op.eval_str(temp) == 1 or data_op.eval_str(temp) == 2 ):
        lmp_dic['temp'] = float(temp)
      else:
        log_info.log_error('Input error: temp should be float or list of float, please check or reset deepff/lammps/temp')
        exit()
  else:
    lmp_dic['temp'] = '300.0'

  if ( 'pres' in lmp_dic.keys() ):
    pres = lmp_dic['pres']
    if ( isinstance(pres, list) ):
      if ( all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in pres)):
        lmp_dic['pres'] = [float(i) for i in pres]
      else:
        log_info.log_error('Input error: multipole pressure should be list of float, please check or reset deepff/lammps/pres')
        exit()
    else:
      if ( data_op.eval_str(pres) == 1 or data_op.eval_str(pres) == 2 ):
        lmp_dic['pres'] = float(pres)
      else:
        log_info.log_error('Input error: pres shoule be float or list of float, please check deepff/lammps/pres')
        exit()
  else:
    lmp_dic['pres'] = '1.0'

  sys_num = 0
  for key in lmp_dic.keys():
    if ( 'system' in key ):
      sys_num = sys_num+1

      if ( 'box' in lmp_dic[key].keys() ):
        box_file = lmp_dic[key]['box']
        if ( os.path.exists(os.path.abspath(box_file)) ):
          lmp_dic[key]['box'] = os.path.abspath(box_file)
        else:
          log_info.log_error('Input error: %s file does not exist' %(box_file))
          exit()

      if ( 'coord' in lmp_dic[key].keys() ):
        coord_file = lmp_dic[key]['coord']
        if ( os.path.exists(os.path.abspath(coord_file)) ):
          lmp_dic[key]['coord'] = os.path.abspath(coord_file)
        else:
          log_info.log_error('Input error: %s file does not exist' %(coord_file))
          exit()

      valid_md_type = ['nve', 'nvt', 'npt']
      if ( 'md_type' in lmp_dic[key].keys() ):
        md_type = lmp_dic[key]['md_type']
        if ( md_type in valid_md_type ):
          pass
        else:
          log_info.log_error('Input error: only nve, nvt and npt are supportted for md_type, please check or reset deepff/lammps/system/md_type')
          exit()
      else:
        lmp_dic[key]['md_type'] = 'nvt'

      if ( 'use_metad' in lmp_dic[key].keys() ):
        use_metad = lmp_dic[key]['use_metad']
        use_metad_bool = data_op.str_to_bool(use_metad)
        if ( isinstance(use_metad_bool, bool) ):
          lmp_dic[key]['use_metad'] = use_metad_bool
        else:
          log_info.log_out('Input error: use_metad should be bool, please check or set deepff/lammps/system/use_metad')
          exit()
      else:
        lmp_dic[key]['use_metad'] = False

      if lmp_dic[key]['use_metad'] :
        if ( 'plumed_file' in lmp_dic[key].keys() ):
          plumed_file = lmp_dic[key]['plumed_file']
          if ( os.path.exists(os.path.abspath(plumed_file)) ):
            lmp_dic[key]['plumed_file'] = os.path.abspath(plumed_file)
          else:
            log_info.log_error('Input error: %s file does not exist, please check deepff/lammps/system/plumed_file' %(plumed_file))
            exit()
        else:
          log_info.log_error('Input error: as user want to use plumed, but no plumed_file, please check deepff/lammps/system/plumed_file')
          exit()

  if ( sys_num == 0 ):
    log_info.log_error('Input error: no system for lammps calculation, please set deepff/lammps/system')
    exit()

  #Check parameters for model_devi
  if ( 'choose_new_data_num_limit' in model_devi_dic.keys() ):
    choose_new_data_num_limit = model_devi_dic['choose_new_data_num_limit']
    if ( data_op.eval_str(choose_new_data_num_limit) == 1 ):
      model_devi_dic['choose_new_data_num_limit'] = int(choose_new_data_num_limit)
    else:
      log_info.log_error('Input error: choose_new_data_num_limit should be integer, please check or reset deepff/model_devi/choose_new_data_num_limit')
      exit()
  else:
    model_devi_dic['choose_new_data_num_limit'] = 100

  if ( 'force_conv' in model_devi_dic.keys() ):
    force_conv = model_devi_dic['force_conv']
    if ( data_op.eval_str(force_conv) == 1 or data_op.eval_str(force_conv) == 2 ):
      model_devi_dic['force_conv'] = float(force_conv)
    else:
      log_info.log_error('Input error: force_conv should be integer or float, please check or set deepff/model_devi/force_conv')
      exit()
  else:
    model_devi_dic['force_conv'] = 0.05

  if ( 'max_iter' in model_devi_dic.keys() ):
    max_iter = model_devi_dic['max_iter']
    if ( data_op.eval_str(max_iter) == 1 ):
      model_devi_dic['max_iter'] = int(max_iter)
    else:
      log_info.log_error('Input error: max_iter should be integer, please check or reset deepff/model_devi/max_iter')
      exit()
  else:
    model_devi_dic['max_iter'] = 100

  if ( 'restart_iter' in model_devi_dic.keys() ):
    restart_iter = model_devi_dic['restart_iter']
    if ( data_op.eval_str(restart_iter) == 1 ):
      model_devi_dic['restart_iter'] = int(restart_iter)
    else:
      log_info.log_error('Input error: restart_iter should be integer, please check or reset deepff/model_devi/restart_iter')
      exit()
  else:
    model_devi_dic['restart_iter'] = 0

  if ( 'restart_data_num' in model_devi_dic.keys() ):
    restart_data_num = model_devi_dic['restart_data_num']
    if ( data_op.eval_str(restart_data_num) == 1 ):
      model_devi_dic['restart_data_num'] = int(restart_data_num)
    else:
      log_info.log_error('Input error: restart_data_num should be integer, please check or reset deepff/model_devi/restart_data_num')
      exit()
  else:
    model_devi_dic['restart_data_num'] = 0

  if ( 'restart_stage' in model_devi_dic.keys() ):
    restart_stage = model_devi_dic['restart_stage']
    if ( data_op.eval_str(restart_stage) == 1 ):
      model_devi_dic['restart_stage'] = int(restart_stage)
    else:
      log_info.log_error('Input error: restart_stage should be integer, please check or reset deepff/model_devi/restart_stage')
      exit()
  else:
    model_devi_dic['restart_stage'] = 0

  #Check parameters for CP2K
  #For multi-system, we need multi cp2k input files.
  cp2k_inp_file_tot = []
  if ( 'cp2k_inp_file' in cp2k_dic.keys() ):
    cp2k_inp_file = cp2k_dic['cp2k_inp_file']
    if ( isinstance(cp2k_inp_file, list) ):
      if ( len(cp2k_inp_file) != sys_num ):
        log_info.log_error('Input error: the number of cp2k input file should be equal to the number of systems, please check or reset deepff/cp2k/cp2k_inp_file')
        exit()
      else:
        for inp_file in cp2k_inp_file:
          if ( os.path.exists(os.path.abspath(inp_file)) ):
            cp2k_inp_file_tot.append(os.path.abspath(inp_file))
          else:
            log_info.log_error('%s file does not exist' %(inp_file))
            exit()
    elif ( isinstance(cp2k_inp_file, str) ):
      if ( os.path.exists(os.path.abspath(cp2k_inp_file)) ):
        cp2k_inp_file_tot.append(os.path.abspath(cp2k_inp_file))
      else:
        log_info.log_error('%s file does not exist' %(cp2k_inp_file))
        exit()
    cp2k_dic['cp2k_inp_file'] = cp2k_inp_file_tot
  else:
    for i in range(sys_num):
      cp2k_inp_file_tot.append('none')
    cp2k_dic['cp2k_inp_file'] = cp2k_inp_file_tot
    if ( 'basis_set_file_name' in cp2k_dic.keys() ):
      basis_set_file_name = cp2k_dic['basis_set_file_name']
      if ( os.path.exists(os.path.abspath(basis_set_file_name)) ):
        cp2k_dic['basis_set_file_name'] = os.path.abspath(basis_set_file_name)
      else:
        log_info.log_error('%s file does not exist' %(basis_set_file_name))
        exit()
    else:
      log_info.log_error('Input error: no basis_set_file_name, please set deepff/cp2k/basis_set_file_name')
      exit()

    if ( 'potential_file_name' in cp2k_dic.keys() ):
      potential_file_name = cp2k_dic['potential_file_name']
      if ( os.path.exists(os.path.abspath(potential_file_name)) ):
        cp2k_dic['potential_file_name'] = os.path.abspath(potential_file_name)
      else:
        log_info.log_error('%s file does not exist' %(potential_file_name))
        exit()
    else:
      log_info.log_error('Input error: no potential_file_name, please set deepff/cp2k/potential_file_name')
      exit()

    if ( 'use_sr_basis' in cp2k_dic.keys() ):
      use_sr_basis = cp2k_dic['use_sr_basis']
      use_sr_basis__bool = data_op.str_to_bool(use_sr_basis)
      if ( isinstance(use_sr_basis_bool, bool) ):
        cp2k_dic['use_sr_basis'] = use_sr_basis_bool
      else:
        log_info.log_error('Input error: use_sr_basis should be bool, please check or reset deepff/cp2k/use_sr_basis')
        exit()
    else:
      cp2k_dic['use_sr_basis'] = False

    if ( 'basis_level' in cp2k_dic.keys() ):
      basis_level = cp2k_dic['basis_level']
      if ( basis_level == 'svp' or basis_level == 'dzvp' or basis_level == 'tzvp' or basis_level == 'tzv2p'):
        pass
      else:
        log_info.log_error('Input error: %s basis set is not surported!' %(basis_level))
        exit()
    else:
      cp2k_dic['basis_level'] = 'dzvp'

    periodic_valid = ['NONE', 'X', 'XY', 'XYZ', 'XZ', 'Y', 'YZ', 'Z']
    if ( 'poisson_periodic' in cp2k_dic.keys() ):
      poisson_periodic = cp2k_dic['poisson_periodic']
      if ( poisson_periodic.upper() in periodic_valid ):
        cp2k_dic['poisson_periodic'] = poisson_periodic.upper()
      else:
        log_info.log_error('Input error: poisson_periodic %s is not supported, please check' %(poisson_periodic))
    else:
      cp2k_dic['poisson_periodic'] = 'XYZ'

    if ( 'cell_periodic' in cp2k_dic.keys() ):
      cell_periodic = cp2k_dic['cell_periodic']
      if ( cell_periodic.upper() in periodic_valid ):
        cp2k_dic['cell_periodic'] = cell_periodic.upper()
      else:
        log_info.log_error('Input error: cell_periodic %s is not supported, please check' %(cell_periodic))
    else:
      cp2k_dic['cell_periodic'] = 'XYZ'

    if ( 'charge' in cp2k_dic.keys() ):
      charge = cp2k_dic['charge']
      if ( data_op.eval_str(charge) == 1 ):
        pass
      else:
        log_info.log_error('Input error: charge should be integer, please check or reset deepff/cp2k/charge')
        exit()
    else:
      log_info.log_error('Input error: no charge, please set deepff/cp2k/charge')
      exit()

    if ( 'multiplicity' in cp2k_dic.keys() ):
      multiplicity = cp2k_dic['multiplicity']
      if ( data_op.eval_str(multiplicity) == 1 ):
        pass
      else:
        log_info.log_error('Input error: multiplicity wrong should be integer, please check or reset deepff/cp2k/multiplicity')
        exit()
    else:
      log_info.log_error('Input error: no multiplicity, please set deepff/cp2k/multiplicity')
      exit()

    if ( 'cutoff' in cp2k_dic.keys() ):
      cutoff = cp2k_dic['cutoff']
      if ( data_op.eval_str(cutoff) == 1 or data_op.eval_str(cutoff) == 2 ):
        pass
      else:
        log_info.log_error('Input error: cutoff should be float or integer, please check deepff/cp2k/cutoff')
        exit()
    else:
      cp2k_dic['cutoff'] = '400'

    functional_lib = ['PBE', 'B3LYP', 'TPSS']
    if ( 'xc_functional' in cp2k_dic.keys() ):
      xc_functional = cp2k_dic['xc_functional']
      if ( xc_functional in functional_lib ):
        pass
      else:
        log_info.log_error('Input error: %s functional is not suported for xc functional' %(xc_functional))
        exit()
    else:
      cp2k_dic['xc_functional'] = 'PBE'

    if ( 'dftd3' in cp2k_dic.keys() ):
      dftd3 = cp2k_dic['dftd3']
      dftd3_bool = data_op.str_to_bool(dftd3)
      if ( isinstance(dftd3_bool, bool) ):
        cp2k_dic['dftd3'] = dftd3_bool
      else:
        log_info.log_error('Input error: dftd3 should be bool, please check or check or reset deepff/cp2k/dftd3')
        exit()
    else:
      cp2k_dic['dftd3'] = False

    if cp2k_dic['dftd3']:
      if ( 'dftd3_file' in cp2k_dic.keys() ):
        dftd3_file = cp2k_dic['dftd3_file']
        if ( os.path.exists(os.path.abspath(dftd3_file)) ):
          cp2k_dic['dftd3_file'] = os.path.abspath(dftd3_file)
        else:
          log_info.log_error('%s file does not exist' %(os.path.abspath(dftd3_file)))
          exit()
      else:
        log_info.log_error('Input error: no dftd3 file, please set deepff/cp2k/dftd3_file')

  if ( 'use_prev_wfn' in cp2k_dic.keys() ):
    use_prev_wfn = cp2k_dic['use_prev_wfn']
    use_prev_wfn_bool = data_op.str_to_bool(use_prev_wfn)
    if ( isinstance(use_prev_wfn_bool, bool) ):
      cp2k_dic['use_prev_wfn'] = use_prev_wfn_bool
    else:
      log_info.log_error('Input error: use_prev_wfn should be bool, please check or reset deepff/cp2k/use_prev_wfn')
      exit()
  else:
    cp2k_dic['use_prev_wfn'] = False

  #Check parameters for environ
  if ( 'cp2k_exe' in environ_dic.keys() ):
    cp2k_exe = environ_dic['cp2k_exe']
    if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
      environ_dic['cp2k_exe'] = os.path.abspath(cp2k_exe)
    else:
      log_info.log_error('Input error: cp2k executable file does not exist, please check or set deepff/environ/cp2k_exe')
      exit()
  else:
    log_info.log_error('Input error: no cp2k executable file, please set deepff/environ/cp2k_exe')
    exit()

  if ( 'cp2k_env_file' in environ_dic.keys() ):
    cp2k_env_file = environ_dic['cp2k_env_file']
    if ( os.path.exists(os.path.abspath(cp2k_env_file)) ):
      environ_dic['cp2k_env_file'] = os.path.abspath(cp2k_env_file)
    else:
      log_info.log_error('Input error: cp2k environment file does not exist, please check or set deepff/environ/cp2k_env_file')
      exit()
  else:
    log_info.log_error('Input error: no cp2k environment file, please set deepff/environ/cp2k_env_file')
    exit()

  if ( 'cuda_dir' in environ_dic.keys() ):
    cuda_dir = environ_dic['cuda_dir']
  else:
    environ_dic['cuda_dir'] = 'none'

  if ( environ_dic['cuda_dir'] != 'none' ):
    if ( os.path.exists(os.path.abspath(cuda_dir)) ):
      environ_dic['cuda_dir'] = os.path.abspath(cuda_dir)
    else:
      log_info.log_error('Input error: cuda directory does not exist, please check or set deepff/environ/cuda_dir')
      exit()

  if ( 'parallel_exe' in environ_dic.keys() ):
    parallel_exe = environ_dic['parallel_exe']
    if ( os.path.exists(os.path.abspath(parallel_exe)) ):
      environ_dic['parallel_exe'] = os.path.abspath(parallel_exe)
    else:
      log_info.log_error('Input error: parallel executable file does not exist, please check or set deepff/environ/parallel_exe')
      exit()
  else:
    log_info.log_error('Input error: no cp2k parallel file, please set deepff/environ/parallel_exe')
    exit()

  if ( 'cp2k_job_per_node' in environ_dic.keys() ):
    cp2k_job_per_node = environ_dic['cp2k_job_per_node']
    if ( data_op.eval_str(cp2k_job_per_node) == 1 ):
      environ_dic['cp2k_job_per_node'] = int(cp2k_job_per_node)
    else:
      log_info.log_error('Input error: cp2k_job_per_node should be integer, please check or reset deepff/environ/cp2k_job_per_node')
      exit()
  else:
    environ_dic['cp2k_job_per_node'] = 1

  if ( 'lmp_mpi_num' in environ_dic.keys() ):
    lmp_mpi_num = environ_dic['lmp_mpi_num']
    if ( data_op.eval_str(lmp_mpi_num) == 1 ):
      environ_dic['lmp_mpi_num'] = int(lmp_mpi_num)
    else:
      log_info.log_error('Input error: lmp_mpi_num should be integer, please check or reset deepff/environ/lmp_mpi_num')
      exit()
  else:
    environ_dic['lmp_mpi_num'] = multiprocessing.cpu_count()

  if ( 'lmp_openmp_num' in environ_dic ):
    lmp_openmp_num = environ_dic['lmp_openmp_num']
    if ( data_op.eval_str(lmp_openmp_num) == 1 ):
      environ_dic['lmp_openmp_num'] = int(lmp_openmp_num)
    else:
      log_info.log_error('Input error: lmp_openmp_num should be integer, please check or reset deepff/environ/lmp_openmp_num')
      exit()
  else:
    environ_dic['lmp_openmp_num'] = 1

  return deepmd_dic, lmp_dic, cp2k_dic, model_devi_dic, environ_dic

def write_restart_inp(inp_file_name, restart_iter, restart_stage, tot_data_num, work_dir):

  '''
  check_inp: write restart input file for deepff

  Args:
    inp_file_name: string
      inp_file_name is the name of cp2k_kit deepff file.
    restart_iter: int
      restart_iter is the iteration number of restart.
    restart_stage: int
      restart_stage is the stage of restart.
    tot_data_num: int
      tot_data_num is the total data number.
    work_dir: string
      work_dir is working directory of CP2K_kit.
  Returns:
    none
  '''

  i=0
  while True:
    inp_tmp_file_name = ''.join((work_dir, '/input_tmp', str(i), '.inp'))
    if ( os.path.exists(inp_tmp_file_name) ):
      i=i+1
    else:
      break

  cmd = "cp %s %s" %(inp_file_name, inp_tmp_file_name)
  call.call_simple_shell(work_dir, cmd)

  restart_inp_file_name = ''.join((work_dir, '/CP2K_KIT.restart'))
  whole_line_num = len(open(inp_tmp_file_name).readlines())

  line_1_num = file_tools.grep_line_num("'&model_devi'", inp_tmp_file_name, work_dir)[0]
  line_2_num = file_tools.grep_line_num("'choose_new_data_num_limit'", inp_tmp_file_name, work_dir)
  line_3_num = file_tools.grep_line_num("'model_devi_freq'", inp_tmp_file_name, work_dir)
  line_4_num = file_tools.grep_line_num("'force_conv'", inp_tmp_file_name, work_dir)
  line_5_num = file_tools.grep_line_num("'max_iter'", inp_tmp_file_name, work_dir)
  line_6_num = file_tools.grep_line_num("'&end model_devi'", inp_tmp_file_name, work_dir)[0]

  restart_inp_file = open(restart_inp_file_name, 'w')
  for i in range(line_1_num):
    line_i = linecache.getline(inp_tmp_file_name, i+1)
    restart_inp_file.write(line_i)

  for i in [line_2_num, line_3_num, line_4_num, line_5_num]:
    if ( i != 0 ):
      line_i = linecache.getline(inp_tmp_file_name, i[0])
      restart_inp_file.write(line_i)

  restart_inp_file.write('    restart_iter %d\n' %(restart_iter))
  restart_inp_file.write('    restart_data_num %d\n' %(tot_data_num))
  restart_inp_file.write('    restart_stage %d\n' %(restart_stage))

  for i in range(line_6_num, whole_line_num+1, 1):
    line_i = linecache.getline(inp_tmp_file_name, i)
    restart_inp_file.write(line_i)

  linecache.clearcache()
  restart_inp_file.close()

  cmd = "rm %s" %(inp_tmp_file_name)
  call.call_simple_shell(work_dir, cmd)
