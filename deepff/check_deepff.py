#! /usr/env/bin python

import os
import copy
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op

def check_inp(deepmd_dic, lmp_dic, cp2k_dic, force_eval_dic, environ_dic, proc_num):

  '''
  check_inp : check the deepff input file

  Args :
    deepmd_dic : dictionary
      deepmd_dic contains keywords used in deepmd.
    lammps_dic : dictionary
      lammpd_dic contains keywords used in lammps.
    cp2k_dic : dictionary
      cp2k_dic contains keywords used in cp2k.
    force_eval_dic : dictionary
      force_eval contains keywords used in force_eval.
    environ_dic : dictionary
      environ_dic contains keywords used in environment.
  Returns :
    new_deepmd_dic : dictionary
      new_deepmd_dic is the revised deepmd_dic.
    new_lmp_dic : dictionary
      new_lmp_dic is the revised lmp_dic.
    new_cp2k_dic : dictionary
      new_cp2k_dic is the revised cp2k_dic.
    new_force_eval_dic : dictionary
      new_force_eval_dic is the revised force_eval_dic.
    new_environ_dic : dictionary
      new_environ_dic is the revised environ_dic.
  '''

  new_deepmd_dic = copy.deepcopy(deepmd_dic)
  new_lmp_dic = copy.deepcopy(lmp_dic)
  new_cp2k_dic = copy.deepcopy(cp2k_dic)
  new_force_eval_dic = copy.deepcopy(force_eval_dic)
  new_environ_dic = copy.deepcopy(environ_dic)

  #Check parameters for deepmd-kit
  if ( 'model' not in new_deepmd_dic.keys() ):
    log_info.log_error('Input error: no model found, please check or set deepff/deepmd/model')
    exit()
  else:
    if ( 'type_map' in new_deepmd_dic['model'].keys() ):
      if ( all(data_op.eval_str(i) == 0 for i in new_deepmd_dic['model']['type_map']) ):
        pass
      else:
        log_info.log_error('Input error: type_map wrong, please check deepff/deepmd/model/type_map')
        exit()
    else:
      log_info.log_error('Input error: no type_map found, please check or set deepff/deepmd/model/type_map')
      exit()
    if ( 'descriptor' not in new_deepmd_dic['model'].keys() ):
      log_info.log_error('Input error: no descriptor found, please check or set deepff/deepmd/model/descriptor')
      exit()
    else:
      if ( 'sel' in new_deepmd_dic['model']['descriptor'].keys() ):
        sel = new_deepmd_dic['model']['descriptor']['sel']
        if ( all(data_op.eval_str(i) == 1 for i in sel) ):
          new_deepmd_dic['model']['descriptor']['sel'] = [int(x) for x in sel]
        else:
          log_info.log_error('Input error: sel wrong, please check deepff/deepmd/model/descriptor/sel')
          exit()
      else:
        log_info.log_error('Input error: no sel, please check or set deepff/deepmd/model/descriptor/sel')
        exit()

      if ( 'rcut_smth' in new_deepmd_dic['model']['descriptor'].keys() ):
        rcut_smth = new_deepmd_dic['model']['descriptor']['rcut_smth']
        if ( data_op.eval_str(rcut_smth) == 1 or data_op.eval_str(rcut_smth) == 2 ):
          new_deepmd_dic['model']['descriptor']['rcut_smth'] = float(rcut_smth)
        else:
          log_info.log_error('Input error: rcut_smth error, please check deepff/deepmd/model/descriptor/rcut_smth')
          exit()
      else:
        new_deepmd_dic['model']['descriptor']['rcut_smth'] = 5.0

      if ( 'rcut' in new_deepmd_dic['model']['descriptor'].keys() ):
        rcut = new_deepmd_dic['model']['descriptor']['rcut']
        if ( data_op.eval_str(rcut) == 1 or data_op.eval_str(rcut) == 2 ):
          new_deepmd_dic['model']['descriptor']['rcut'] = float(rcut)
        else:
          log_info.log_error('Input error: rcut error, please check deepff/deepmd/model/descriptor/rcut')
          exit()
      else:
        new_deepmd_dic['model']['descriptor']['rcut'] = 6.0

      if ( 'neuron' in new_deepmd_dic['model']['descriptor'].keys() ):
        neuron_encode = new_deepmd_dic['model']['descriptor']['neuron']
        if ( all(data_op.eval_str(i) == 1 for i in neuron_encode) ):
          new_deepmd_dic['model']['descriptor']['neuron'] = [int(x) for x in neuron_encode]
        else:
          log_info.log_error('Input error: neuron error, please check deepff/deepmd/model/descriptor/neuron')
          exit()
      else:
        new_deepmd_dic['model']['descriptor']['neuron'] = [25, 50, 100]

      if ( 'resnet_dt' in new_deepmd_dic['model']['descriptor'].keys() ):
        resnet_dt = new_deepmd_dic['model']['descriptor']['resnet_dt']
        resnet_dt_bool = data_op.str_to_bool(resnet_dt)
        if ( isinstance(resnet_dt_bool, bool) ):
          new_deepmd_dic['model']['descriptor']['resnet_dt'] = resnet_dt_bool
        else:
          log_info.log_error('Input error: resnet_dt error, please check deepff/deepmd/model/descriptor/resnet_dt')
          exit()
      else:
        new_deepmd_dic['model']['descriptor']['resnet_dt'] = False

      if ( 'axis_neuron' in new_deepmd_dic['model']['descriptor'].keys() ):
        axis_neuron =new_deepmd_dic['model']['descriptor']['axis_neuron']
        if ( data_op.eval_str(axis_neuron) == 1 ):
          new_deepmd_dic['model']['descriptor']['axis_neuron'] = int(axis_neuron)
        else:
          log_info.log_error('Input error: axis_neuron error, please check deepff/deepmd/model/descriptor/axis_neuron')
          exit()
      else:
        new_deepmd_dic['model']['descriptor']['axis_neuron'] = 16

    if ( 'fitting_net' not in new_deepmd_dic['model'].keys() ):
      log_info.log_error('Input error: no fitting_net found, please check or set deepff/deepmd/model/fitting_net')
      exit()
    else:
      if ( 'resnet_dt' in new_deepmd_dic['model']['fitting_net'].keys() ):
        resnet_dt = new_deepmd_dic['model']['fitting_net']['resnet_dt']
        resnet_dt_bool = data_op.str_to_bool(resnet_dt)
        if ( isinstance(resnet_dt_bool, bool) ):
          new_deepmd_dic['model']['fitting_net']['resnet_dt'] = resnet_dt_bool
        else:
          log_info.log_error('Input error: resnet_dt error, please check deepff/deepmd/model/fitting_net/resnet_dt')
          exit()
      else:
        new_deepmd_dic['model']['fitting_net']['resnet_dt'] = True

  if ( 'learning_rate' not in new_deepmd_dic.keys() ):
    log_info.log_error('Input error: no learning_rate found, please check or set deepff/deepmd/learning_rate')
    exit()
  else:
    if ( 'type' in new_deepmd_dic['learning_rate'].keys() ):
      decay_type = new_deepmd_dic['learning_rate']['type']
      if ( data_op.eval_str(decay_type) == 0 ):
        pass
      else:
        log_info.log_error('Input error: type error, please check deepff/deepmd/learning_rate/type')
        eixt()
    else:
      new_deepmd_dic['learning_rate']['type'] = 'exp'

    if ( 'decay_steps' in new_deepmd_dic['learning_rate'].keys() ):
      decay_steps = new_deepmd_dic['learning_rate']['decay_steps']
      if ( data_op.eval_str(decay_steps) == 1 ):
        new_deepmd_dic['learning_rate']['decay_steps'] = int(decay_steps)
      else:
        log_info.log_error('Input error: decay_steps error, please check deepff/deepmd/learning_rate/decay_steps')
        eixt()
    else:
      new_deepmd_dic['learning_rate']['decay_steps'] = 5000

    if ( 'start_lr' in new_deepmd_dic['learning_rate'].keys() ):
      start_lr = new_deepmd_dic['learning_rate']['start_lr']
      if ( data_op.eval_str(start_lr) == 2 ):
        new_deepmd_dic['learning_rate']['start_lr'] = float(start_lr)
      else:
        log_info.log_error('Input error: start_lr error, please check deepff/deepmd/learning_rate/start_lr')
        exit()
    else:
      new_deepmd_dic['learning_rate']['start_lr'] = 0.001

    if ( 'stop_lr' in new_deepmd_dic['learning_rate'].keys() ):
      stop_lr = new_deepmd_dic['learning_rate']['stop_lr']
      if ( data_op.eval_str(stop_lr) == 2 ):
        new_deepmd_dic['learning_rate']['stop_lr'] = float(stop_lr)
      else:
        log_info.log_error('Input error: stop_lr error, please check deepff/deepmd/learning_rate/stop_lr')
        exit()
    else:
      new_deepmd_dic['learning_rate']['stop_lr'] = 1e-8

  if ( 'loss' not in new_deepmd_dic.keys() ):
    log_info.log_error('Input error: no loss found, please check or set deepff/deepmd/loss')
    exit()
  else:
    if ( 'start_pref_e' in new_deepmd_dic['loss'].keys() ):
      start_pref_e = new_deepmd_dic['loss']['start_pref_e']
      if ( data_op.eval_str(start_pref_e) == 1 or data_op.eval_str(start_pref_e) == 2 ):
        new_deepmd_dic['loss']['start_pref_e'] = float(start_pref_e)
      else:
        log_info.log_error('Input error: start_pref_e error, please check deepff/deepmd/loss/start_pref_e')
        exit()
    else:
      new_deepmd_dic['loss']['start_pref_e'] = 0.02

    if ( 'limit_pref_e' in new_deepmd_dic['loss'].keys() ):
      limit_pref_e = new_deepmd_dic['loss']['limit_pref_e']
      if ( data_op.eval_str(limit_pref_e) == 1 or data_op.eval_str(limit_pref_e) == 2 ):
        new_deepmd_dic['loss']['limit_pref_e'] = float(limit_pref_e)
      else:
        log_info.log_error('Input error: limit_pref_e error, please check deepff/deepmd/loss/limit_pref_e')
        exit()
    else:
      new_deepmd_dic['loss']['limit_pref_e'] = 1.0

    if ( 'start_pref_f' in new_deepmd_dic['loss'].keys() ):
      start_pref_f = new_deepmd_dic['loss']['start_pref_f']
      if ( data_op.eval_str(start_pref_f) == 1 or data_op.eval_str(start_pref_f) == 2 ):
        new_deepmd_dic['loss']['start_pref_f'] = float(start_pref_f)
      else:
        log_info.log_error('Input error: start_pref_f error, please check deepff/deepmd/loss/start_pref_f')
        exit()
    else:
      new_deepmd_dic['loss']['start_pref_f'] = 1000.0

    if ( 'limit_pref_f' in new_deepmd_dic['loss'].keys() ):
      limit_pref_f = new_deepmd_dic['loss']['limit_pref_f']
      if ( data_op.eval_str(limit_pref_f) == 1 or data_op.eval_str(limit_pref_f) == 2 ):
        new_deepmd_dic['loss']['limit_pref_f'] = float(limit_pref_f)
      else:
        log_info.log_error('Input error: limit_pref_f error, please check deepff/deepmd/loss/limit_pref_f')
        exit()
    else:
      new_deepmd_dic['loss']['limit_pref_f'] = 1.0

    if ( 'start_pref_v' in new_deepmd_dic['loss'].keys() ):
      start_pref_v = new_deepmd_dic['loss']['start_pref_v']
      if ( data_op.eval_str(start_pref_v) == 1 or data_op.eval_str(start_pref_v) == 2 ):
        new_deepmd_dic['loss']['start_pref_v'] = float(start_pref_v)
      else:
        log_info.log_error('Input error: start_pref_v error, please check deepff/deepmd/loss/start_pref_v')
        exit()
    else:
      new_deepmd_dic['loss']['start_pref_v'] = 0.0

    if ( 'limit_pref_v' in new_deepmd_dic['loss'].keys() ):
      limit_pref_v = new_deepmd_dic['loss']['limit_pref_v']
      if ( data_op.eval_str(limit_pref_v) == 1 or data_op.eval_str(limit_pref_v) == 2 ):
        new_deepmd_dic['loss']['limit_pref_v'] = float(limit_pref_v)
      else:
        log_info.log_error('Input error: limit_pref_v error, please check deepff/deepmd/loss/limit_pref_v')
        exit()
    else:
      new_deepmd_dic['loss']['limit_pref_v'] = 0.0

  if ( 'training' not in new_deepmd_dic.keys() ):
    log_info.log_error('Input error: no training found, please check or set deepff/deepmd/training')
    exit()
  else:
    for i in new_deepmd_dic['training'].keys():
      if ( 'system' in i ):
        if ( 'directory' in new_deepmd_dic['training'][i] ):
          directory = new_deepmd_dic['training'][i]['directory']
          if ( os.path.exists(os.path.abspath(directory)) ):
            pass
          else:
            log_info.log_error('Input error: %s does not exist, please check deepff/deepmd/training/system/directory' %(directory))
            eixt()
        else:
          log_info.log_error('Input error: no directory, please check deepff/deepmd/training/system/directory')
          exit()
        if ( 'proj_name' in new_deepmd_dic['training'][i]):
          pass
        else:
          log_info.log_error('Input error: no proj_name, please check deepff/deepmd/training/system/proj_name')
          exit()
        if ( 'start_frame' in new_deepmd_dic['training'][i] ):
          start_frame = new_deepmd_dic['training'][i]['start_frame']
          if ( data_op.eval_str(start_frame) == 1 ):
            new_deepmd_dic['training'][i]['start_frame'] = int(start_frame)
          else:
            log_info.log_error('Input error: start_frame error, please check deepff/deepmd/training/system/start_frame')
            exit()
        else:
          new_deepmd_dic['training'][i]['start_frame'] = 0
        if ( 'end_frame' in new_deepmd_dic['training'][i] ):
          end_frame = new_deepmd_dic['training'][i]['end_frame']
          if ( data_op.eval_str(end_frame) == 1 ):
            new_deepmd_dic['training'][i]['end_frame'] = int(end_frame)
          else:
            log_info.log_error('Input error: end_frame error, please check deepff/deepmd/training/system/end_frame')
            exit()
        else:
          new_deepmd_dic['training'][i]['end_frame'] = 0
        if ( 'choosed_frame_num' in new_deepmd_dic['training'][i] ):
          choosed_frame_num = new_deepmd_dic['training'][i]['choosed_frame_num']
          if ( data_op.eval_str(choosed_frame_num) == 1 ):
            new_deepmd_dic['training'][i]['choosed_frame_num'] = int(choosed_frame_num)
          else:
            log_info.log_error('Input error: choosed_frame_num error, please check deepff/deepmd/training/system/choosed_frame_num')
            exit()
        else:
          new_deepmd_dic['training'][i]['choosed_frame_num'] = 0
        if ( 'set_parts' in new_deepmd_dic['training'][i] ):
          set_parts = new_deepmd_dic['training'][i]['set_parts']
          if ( data_op.eval_str(set_parts) == 1 ):
            new_deepmd_dic['training'][i]['set_parts'] = int(set_parts)
          else:
            log_info.log_error('Input error: set_parts error, please check deepff/deepmd/training/system/set_parts')
            exit()
        else:
          new_deepmd_dic['training'][i]['set_parts'] = 1

    if ( 'set_data_dir' in new_deepmd_dic['training'].keys() ):
      set_data_dir = new_deepmd_dic['training']['set_data_dir']
      if ( os.path.exists(os.path.abspath(set_data_dir)) ):
        pass
      else:
        log_info.log_error('Input error: %s directory does not exist, please check or set deepff/deepmd/training/set_data_dir')
        exit()

    if ( 'model_type' in new_deepmd_dic['training'].keys() ):
      model_type = new_deepmd_dic['training']['model_type']
      if ( model_type == 'use_seed' or model_type == 'use_node' ):
        pass
      else:
        log_info.log_error('Input error: model_type should be use_seed or use_node, please check deepff/deepmd/training/model_type')
        exit()
    else:
      new_deepmd_dic['training']['model_type'] = 'use_seed'

    if ( new_deepmd_dic['training']['model_type'] == 'use_seed' ):
      if ( 'seed_num' in new_deepmd_dic['training'].keys() ):
        seed_num = new_deepmd_dic['training']['seed_num']
        if ( data_op.eval_str(seed_num) == 1 ):
          new_deepmd_dic['training']['seed_num'] = int(seed_num)
        else:
          log_info.log_error('Input error: seed_num error, please check deepff/deepmd/training/seed_num')
          exit()
      else:
        new_deepmd_dic['training']['seed_num'] = 2

    if ( 'neuron' in new_deepmd_dic['training'].keys() ):
      neuron_list = new_deepmd_dic['training']['neuron']
      if ( new_deepmd_dic['training']['model_type'] == 'use_node' ):
        neuron = []
        tmp_str = data_op.comb_list_2_str(neuron_list, ' ')
        tmp_list = data_op.str_split(tmp_str, '...')

        for i in range(len(tmp_list)):
          neuron_i = data_op.str_split(tmp_list[i], ' ')
          if ( all(data_op.eval_str(j) == 1 for j in neuron_i) ):
            neuron.append([int(x) for x in neuron_i])
          else:
            log_info.log_error('Input error: neuron error, please check deepff/deepmd/training/neuron')
            exit()
      elif ( new_deepmd_dic['training']['model_type'] == 'use_seed' ):
        if ( all(data_op.eval_str(j) == 1 for j in neuron_list) ):
          neuron = [int(x) for x in neuron_str]
        else:
          log_info.log_error('Input error: neuron error, please check deepff/deepmd/training/neuron')
          exit()
      new_deepmd_dic['training']['neuron'] = neuron
    else:
      log_info.log_error('Input error: no neuron, please check deepff/deepmd/training/neuron')
      exit()

    if ( 'set_prefix' in new_deepmd_dic['training'].keys() ):
      set_prefix = new_deepmd_dic['training']['set_prefix']
      if ( data_op.eval_str(set_prefix) == 0 ):
        pass
      else:
        log_info.log_error('Input error: set_prefix error, please check deepff/deepmd/training/set_prefix')
        exit()
    else:
      log_info.log_error('Input error: no set_prefix, please check or set deepff/deepmd/training/set_prefix')
      exit()

    if ( 'stop_batch' in new_deepmd_dic['training'].keys() ):
      stop_batch = new_deepmd_dic['training']['stop_batch']
      if ( data_op.eval_str(stop_batch) == 1 ):
        new_deepmd_dic['training']['stop_batch'] = int(stop_batch)
      else:
        log_info.log_error('Input error: stop_batch error, please check deepff/deepmd/training/stop_batch')
        exit()
    else:
      new_deepmd_dic['training']['stop_batch'] = 100000

    if ( 'batch_size' in new_deepmd_dic['training'].keys() ):
      batch_size = new_deepmd_dic['training']['batch_size']
      if ( data_op.eval_str(batch_size) == 1 ):
        new_deepmd_dic['training']['batch_size'] = int(batch_size)
      elif ( batch_size == 'auto' ):
        pass
      else:
        log_info.log_error('Input error: batch_size error, please check deepff/deepmd/training/batch_size')
        exit()
    else:
      new_deepmd_dic['training']['batch_size'] = 1

    if ( 'disp_freq' in new_deepmd_dic['training'].keys() ):
      disp_freq = new_deepmd_dic['training']['disp_freq']
      if ( data_op.eval_str(disp_freq) == 1 ):
        new_deepmd_dic['training']['disp_freq'] = int(disp_freq)
      else:
        log_info.log_error('Input error: disp_freq error, please check deepff/deepmd/training/disp_freq')
        exit()
    else:
      new_deepmd_dic['training']['disp_freq'] = 100

    if ( 'numb_test' in new_deepmd_dic['training'].keys() ):
      numb_test = new_deepmd_dic['training']['numb_test']
      if ( data_op.eval_str(numb_test) == 1 ):
        new_deepmd_dic['training']['numb_test'] = int(numb_test)
      else:
        log_info.log_error('Input error: numb_test error, please check deepff/deepmd/training/numb_test')
        exit()
    else:
      new_deepmd_dic['training']['numb_test'] = 10

    if ( 'save_freq' in new_deepmd_dic['training'].keys() ):
      save_freq = new_deepmd_dic['training']['save_freq']
      if ( data_op.eval_str(save_freq) == 1 ):
        new_deepmd_dic['training']['save_freq'] = int(save_freq)
      else:
        log_info.log_error('Input error: save_freq error, please check deepff/deepmd/training/save_freq')
        exit()
    else:
      new_deepmd_dic['training']['save_freq'] = 1000

    if ( 'disp_training' in new_deepmd_dic['training'].keys() ):
      disp_training = new_deepmd_dic['training']['disp_training']
      disp_training_bool = data_op.str_to_bool(disp_training)
      if ( isinstance(disp_training_bool, bool) ):
        new_deepmd_dic['training']['disp_training'] = disp_training_bool
      else:
        log_info.log_error('Input error: disp_training error, please check deepff/deepmd/training/disp_training')
        exit()
    else:
      new_deepmd_dic['training']['disp_training'] = True

    if ( 'disp_file' in new_deepmd_dic['training'].keys() ):
      disp_file = new_deepmd_dic['training']['disp_file']
      if ( data_op.eval_str(disp_file) == 0 ):
        pass
      else:
        log_info.log_error('Input error: disp_file error, please check deepff/deepmd/training/disp_file')
        exit()
    else:
      new_deepmd_dic['training']['disp_file'] = 'lcurve.out'

    if ( 'save_ckpt' in new_deepmd_dic['training'].keys() ):
      save_ckpt = new_deepmd_dic['training']['save_ckpt']
      if ( data_op.eval_str(save_ckpt) == 0 ):
        pass
      else:
        log_info.log_error('Input error: save_ckpt error, please check deepff/deepmd/training/save_ckpt')
        exit()
    else:
      new_deepmd_dic['training']['load_ckpt'] = 'model.ckpt'

    if ( 'load_ckpt' in new_deepmd_dic['training'].keys() ):
      load_ckpt = new_deepmd_dic['training']['load_ckpt']
      if ( data_op.eval_str(load_ckpt) == 0 ):
        pass
      else:
        log_info.log_error('Input error: load_ckpt error, please check deepff/deepmd/training/load_ckpt')
        exit()
    else:
      new_deepmd_dic['training']['load_ckpt'] = 'model.ckpt'

    if ( 'time_training' in new_deepmd_dic['training'].keys() ):
      time_training = new_deepmd_dic['training']['time_training']
      time_training_bool = data_op.str_to_bool(time_training)
      if ( isinstance(time_training_bool, bool) ):
        new_deepmd_dic['training']['time_training'] = time_training_bool
      else:
        log_info.log_error('Input error: time_training error, please check deepff/deepmd/training/time_training')
        exit()
    else:
      new_deepmd_dic['training']['time_training'] = True

    if ( 'profiling' in new_deepmd_dic['training'].keys() ):
      profiling = new_deepmd_dic['training']['profiling']
      profiling_bool = data_op.str_to_bool(profiling)
      if ( isinstance(profiling_bool, bool) ):
        new_deepmd_dic['training']['profiling'] = profiling_bool
      else:
        log_info.log_error('Input error: profiling error, please check deepff/deepmd/training/profiling')
        exit()
    else:
      new_deepmd_dic['training']['profiling'] = False

  #Check parameters for lammps
  if ( 'nsteps' in new_lmp_dic.keys() ):
    nsteps = new_lmp_dic['nsteps']
    if ( data_op.eval_str(nsteps) == 1 ):
      pass
    else:
      log_info.log_error('Input error: nsteps error, please check deepff/lammps/nsteps')
      exit()
  else:
    new_lmp_dic['nsteps'] = '100000'

  if ( 'thermo_freq' in new_lmp_dic.keys() ):
    thermo_freq = new_lmp_dic['thermo_freq']
    if ( data_op.eval_str(thermo_freq) == 1 ):
      pass
    else:
      log_info.log_error('Input error: thermo_freq error, please check deepff/lammps/thermo_freq')
      exit()
  else:
    new_lmp_dic['thermo_freq'] = '10'

  if ( 'dump_freq' in new_lmp_dic.keys() ):
    dump_freq = new_lmp_dic['dump_freq']
    if ( data_op.eval_str(dump_freq) == 1 ):
      pass
    else:
      log_info.log_error('Input error: dump_freq error, please check deepff/lammps/dump_freq')
      exit()
  else:
    new_lmp_dic['dump_freq'] = '10'

  if ( 'time_step' in new_lmp_dic.keys() ):
    time_step = new_lmp_dic['time_step']
    if ( data_op.eval_str(time_step) == 2 ):
      pass
    else:
      log_info.log_error('Input error: time_step error, please check deepff/lammps/time_step')
  else:
    new_lmp_dic['time_step'] = '0.0005'

  if ( 'tau_t' in new_lmp_dic.keys() ):
    tau_t = new_lmp_dic['tau_t']
    if ( data_op.eval_str(tau_t) == 2 ):
      pass
    else:
      log_info.log_error('Input error: tau_t error, please check deepff/lammps/tau_t')
      exit()
  else:
    new_lmp_dic['tau_t'] = '%f' %(float(new_lmp_dic['time_step'])*200)

  if ( 'tau_p' in new_lmp_dic.keys() ):
    tau_p = new_lmp_dic['tau_p']
    if ( data_op.eval_str(tau_p) == 2 ):
      pass
    else:
      log_info.log_error('Input error: tau_p error, please check deepff/lammps/tau_p')
      exit()
  else:
    new_lmp_dic['tau_p'] = '%f' %(float(new_lmp_dic['time_step'])*200)

  if ( 'md_type' in new_lmp_dic.keys() ):
    md_type = new_lmp_dic['md_type']
    if ( md_type == 'nvt' or md_type == 'npt' ):
      pass
    else:
      log_info.log_error('Input error: md_type wrong, please check deepff/lammps/md_type')
      exit()
  else:
    new_lmp_dic['md_type'] = 'nvt'

  if ( 'temp' in new_lmp_dic.keys() ):
    temp = new_lmp_dic['temp']
    if ( isinstance(temp, list) ):
      if ( all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in temp)):
        pass
      else:
        log_info.log_error('Input error: temp wrong, please check deepff/lammps/temp')
        exit()
    else:
      if ( data_op.eval_str(temp) == 1 or data_op.eval_str(temp) == 2 ):
        pass
      else:
        log_info.log_error('Input error: temp wrong, please check deepff/lammps/temp')
        exit()
  else:
    new_lmp_dic['temp'] = '300.0'

  if ( 'pres' in new_lmp_dic.keys() ):
    pres = new_lmp_dic['pres']
    if ( isinstance(pres, list) ):
      if ( all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in pres)):
        pass
      else:
        log_info.log_error('Input error: pres wrong, please check deepff/lammps/pres')
        exit()
    else:
      if ( data_op.eval_str(pres) == 1 or data_op.eval_str(pres) == 2 ):
        pass
      else:
        log_info.log_error('Input error: pres wrong, please check deepff/lammps/pres')
        exit()
  else:
    new_lmp_dic['pres'] = '1.0'

  sys_num = 0
  for key in new_lmp_dic.keys():
    if ( 'system' in key ):
      sys_num = sys_num+1
      if ( 'box' in new_lmp_dic[key] ):
        box_file = new_lmp_dic[key]['box']
        if ( os.path.exists(os.path.abspath(box_file)) ):
          new_lmp_dic[key]['box'] = os.path.abspath(box_file)
        else:
          log_info.log_error('%s file does not exist' %(box_file))
          exit()
      if ( 'coord' in new_lmp_dic[key] ):
        coord_file = new_lmp_dic[key]['coord']
        if ( os.path.exists(os.path.abspath(coord_file)) ):
          new_lmp_dic[key]['coord'] = os.path.abspath(coord_file)
        else:
          log_info.log_error('%s file does not exist' %(coord_file))
          exit()

  if ( sys_num == 0 ):
    log_info.log_error('Input error: no system for lammps calculation, please set deepff/lammps/system')

  #Check parameters for force_eval
  if ( 'choose_new_data_num_limit' in new_force_eval_dic.keys() ):
    choose_new_data_num_limit = new_force_eval_dic['choose_new_data_num_limit']
    if ( data_op.eval_str(choose_new_data_num_limit) == 1 ):
      new_force_eval_dic['choose_new_data_num_limit'] = int(choose_new_data_num_limit)
    else:
      log_info.log_error('Input error: choose_new_data_num_limit wrong, please check or set deepff/force_eval/choose_new_data_num_limit')
      exit()
  else:
    new_force_eval_dic['choose_new_data_num_limit'] = 100

  if ( 'conv_new_data_num' in new_force_eval_dic.keys() ):
    conv_new_data_num = new_force_eval_dic['conv_new_data_num']
    if ( data_op.eval_str(conv_new_data_num) == 1 ):
      new_force_eval_dic['conv_new_data_num'] = int(conv_new_data_num)
    else:
      log_info.log_error('Input error: conv_new_data_num wrong, please check or set deepff/force_eval/conv_new_data_num')
      exit()
  else:
    new_force_eval_dic['conv_new_data_num'] = 5

  if ( 'force_conv' in new_force_eval_dic.keys() ):
    force_conv = new_force_eval_dic['force_conv']
    if ( data_op.eval_str(conv_new_data_num) == 1 or data_op.eval_str(conv_new_data_num) == 2 ):
      new_force_eval_dic['force_conv'] = float(force_conv)
    else:
      log_info.log_error('Input error: force_conv wrong, please check or set deepff/force_eval/force_conv')
      exit()
  else:
    new_force_eval_dic['force_conv'] = 0.05

  if ( 'max_iter' in new_force_eval_dic.keys() ):
    max_iter = new_force_eval_dic['max_iter']
    if ( data_op.eval_str(max_iter) == 1 ):
      new_force_eval_dic['max_iter'] = int(max_iter)
    else:
      log_info.log_error('Input error: max_iter wrong, please check or set deepff/force_eval/max_iter')
      exit()
  else:
    new_force_eval_dic['max_iter'] = 100

  if ( 'restart_iter' in new_force_eval_dic.keys() ):
    restart_iter = new_force_eval_dic['restart_iter']
    if ( data_op.eval_str(restart_iter) == 1 ):
      new_force_eval_dic['restart_iter'] = int(restart_iter)
    else:
      log_info.log_error('Input error: restart_iter wrong, please check or set deepff/force_eval/restart_iter')
      exit()
  else:
    new_force_eval_dic['restart_iter'] = 0

  #Check parameters for CP2K
  if ( 'cp2k_inp_file' in new_cp2k_dic.keys() ):
    if ( len(new_cp2k_dic.keys()) > 1):
      log_info.log_error('Warning: as cp2k_inp_file is set, other keywords are ignored!', 'Warning')
    else:
      cp2k_inp_file = new_cp2k_dic['cp2k_inp_file']
      if ( os.path.exists(os.path.abspath(cp2k_inp_file)) ):
        new_cp2k_dic['cp2k_inp_file'] = os.path.abspath(cp2k_inp_file)
      else:
        log_info.log_error('%s file does not exist' %(cp2k_inp_file))
        exit()
  else:
    new_cp2k_dic['cp2k_inp_file'] = 'none'
    if ( 'basis_set_file_name' in new_cp2k_dic.keys() ):
      basis_set_file_name = new_cp2k_dic['basis_set_file_name']
      if ( os.path.exists(os.path.abspath(basis_set_file_name)) ):
        new_cp2k_dic['basis_set_file_name'] = os.path.abspath(basis_set_file_name)
      else:
        log_info.log_error('%s file does not exist' %(basis_set_file_name))
        exit()
    else:
      log_info.log_error('Input error: no basis_set_file_name, please set deepff/cp2k/basis_set_file_name')
      exit()

    if ( 'potential_file_name' in new_cp2k_dic.keys() ):
      potential_file_name = new_cp2k_dic['potential_file_name']
      if ( os.path.exists(os.path.abspath(potential_file_name)) ):
        new_cp2k_dic['potential_file_name'] = os.path.abspath(potential_file_name)
      else:
        log_info.log_error('%s file does not exist' %(potential_file_name))
        exit()
    else:
      log_info.log_error('Input error: no potential_file_name, please set deepff/cp2k/potential_file_name')
      exit()

    if ( 'basis_level' in new_cp2k_dic.keys() ):
      basis_level = new_cp2k_dic['basis_level']
      if ( basis_level == 'svp' or basis_level == 'dzvp' or basis_level == 'tzvp' or basis_level == 'tzv2p'):
        pass
      else:
        log_info.log_error('Input error: %s basis set is not surported!' %(basis_level))
        exit()
    else:
      new_cp2k_dic['basis_level'] = 'dzvp'

    periodic_valid = ['NONE', 'X', 'XY', 'XYZ', 'XZ', 'Y', 'YZ', 'Z']
    if ( 'poisson_periodic' in new_cp2k_dic.keys() ):
      poisson_periodic = new_cp2k_dic['poisson_periodic']
      if ( poisson_periodic.upper() in periodic_valid ):
        new_cp2k_dic['poisson_periodic'] = poisson_periodic.upper()
      else:
        log_info.log_error('Input error: poisson_periodic %s is not supported, please check' %(poisson_periodic))
    else:
      new_cp2k_dic['poisson_periodic'] = 'XYZ'

    if ( 'cell_periodic' in new_cp2k_dic.keys() ):
      cell_periodic = new_cp2k_dic['cell_periodic']
      if ( cell_periodic.upper() in periodic_valid ):
        new_cp2k_dic['cell_periodic'] = cell_periodic.upper()
      else:
        log_info.log_error('Input error: cell_periodic %s is not supported, please check' %(cell_periodic))
    else:
      new_cp2k_dic['cell_periodic'] = 'XYZ'

    if ( 'charge' in new_cp2k_dic.keys() ):
      charge = new_cp2k_dic['charge']
      if ( data_op.eval_str(charge) == 1 ):
        pass
      else:
        log_info.log_error('Input error: charge wrong, please check deepff/cp2k/charge')
        exit()
    else:
      log_info.log_error('Input error: no charge, please set deepff/cp2k/charge')
      exit()

    if ( 'multiplicity' in new_cp2k_dic.keys() ):
      multiplicity = new_cp2k_dic['multiplicity']
      if ( data_op.eval_str(multiplicity) == 1 ):
        pass
      else:
        log_info.log_error('Input error: multiplicity wrong, please check deepff/cp2k/multiplicity')
        exit()
    else:
      log_info.log_error('Input error: no multiplicity, please set deepff/cp2k/multiplicity')
      exit()

    if ( 'cutoff' in new_cp2k_dic.keys() ):
      cutoff = new_cp2k_dic['cutoff']
      if ( data_op.eval_str(cutoff) == 1 or data_op.eval_str(cutoff) == 2 ):
        pass
      else:
        log_info.log_error('Input error: cutoff wrong, please check deepff/cp2k/cutoff')
        exit()
    else:
      new_cp2k_dic['cutoff'] = '400'

    functional_lib = ['PBE', 'B3LYP', 'TPSS']
    if ( 'xc_functional' in new_cp2k_dic.keys() ):
      xc_functional = new_cp2k_dic['xc_functional']
      if ( xc_functional in function_lib ):
        pass
      else:
        log_info.log_error('Input error: %s functional is not suported' %(xc_functional))
        exit()
    else:
      new_cp2k_dic['xc_functional'] = 'PBE'

    if ( 'dftd3' in new_cp2k_dic.keys() ):
      dftd3 = new_cp2k_dic['dftd3']
      dftd3_bool = data_op.str_to_bool(dftd3)
      if ( isinstance(dftd3_bool, bool) ):
        new_cp2k_dic['dftd3'] = dftd3_bool
      else:
        log_info.log_error('Input error: dftd3 wrong, please check or set deepff/cp2k/dftd3')
        exit()
    else:
      new_cp2k_dic['dftd3'] = False

    if new_cp2k_dic['dftd3']:
      if ( 'dftd3_file' in new_cp2k_dic.keys() ):
        dftd3_file = new_cp2k_dic['dftd3_file']
        if ( os.path.exists(os.path.abspath(dftd3_file)) ):
          new_cp2k_dic['dftd3_file'] = os.path.abspath(dftd3_file)
        else:
          log_info.log_error('%s file does not exist' %(os.path.abspath(dftd3_file)))
          exit()
      else:
        log_info.log_error('Input error: no dftd3 file, please set deepff/cp2k/dftd3_file')

  #Check parameters for environ
  if ( 'cp2k_exe' in new_environ_dic.keys() ):
    cp2k_exe = new_environ_dic['cp2k_exe']
    if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
      new_environ_dic['cp2k_exe'] = os.path.abspath(cp2k_exe)
    else:
      log_info.log_error('Input error: cp2k executable file does not exist, please check or set deepff/environ/cp2k_exe')
      exit()
  else:
    log_info.log_error('Input error: no cp2k executable file, please set deepff/environ/cp2k_exe')
    exit()

  if ( 'cp2k_env_file' in new_environ_dic.keys() ):
    cp2k_env_file = new_environ_dic['cp2k_env_file']
    if ( os.path.exists(os.path.abspath(cp2k_env_file)) ):
      new_environ_dic['cp2k_env_file'] = os.path.abspath(cp2k_env_file)
    else:
      log_info.log_error('Input error: cp2k environment file does not exist, please check or set deepff/environ/cp2k_env_file')
      exit()
  else:
    log_info.log_error('Input error: no cp2k environment file, please set deepff/environ/cp2k_env_file')
    exit()

  if ( 'cuda_dir' in new_environ_dic.keys() ):
    cuda_dir = new_environ_dic['cuda_dir']
  else:
    new_environ_dic['cuda_dir'] = 'none'

  if ( new_environ_dic['cuda_dir'] != 'none' ):
    if ( os.path.exists(os.path.abspath(cuda_dir)) ):
      new_environ_dic['cuda_dir'] = os.path.abspath(cuda_dir)
    else:
      log_info.log_error('Input error: cuda directory does not exist, please check or set deepff/environ/cuda_dir')
      exit()

  if ( 'parallel_exe' in new_environ_dic.keys() ):
    parallel_exe = new_environ_dic['parallel_exe']
    if ( os.path.exists(os.path.abspath(parallel_exe)) ):
      new_environ_dic['parallel_exe'] = os.path.abspath(parallel_exe)
    else:
      log_info.log_error('Input error: parallel executable file does not exist, please check or set deepff/environ/parallel_exe')
      exit()
  else:
    log_info.log_error('Input error: no cp2k parallel file, please set deepff/environ/parallel_exe')
    exit()

  if ( 'cp2k_job_num' in new_environ_dic.keys() ):
    cp2k_job_num = new_environ_dic['cp2k_job_num']
    if ( data_op.eval_str(cp2k_job_num) == 1 ):
      new_environ_dic['cp2k_job_num'] = int(cp2k_job_num)
    else:
      log_info.log_error('Input error: cp2k_job_num wrong, please check or set deepff/environ/cp2k_job_num')
      exit()
  else:
    new_environ_dic['cp2k_job_num'] = 1

  if ( 'lmp_mpi_num' in new_environ_dic.keys() ):
    lmp_mpi_num = new_environ_dic['lmp_mpi_num']
    if ( data_op.eval_str(lmp_mpi_num) == 1 ):
      new_environ_dic['lmp_mpi_num'] = int(lmp_mpi_num)
    else:
      log_info.log_error('Input error: lmp_mpi_num wrong, please check or set deepff/environ/lmp_mpi_num')
      exit()
  else:
    new_environ_dic['lmp_mpi_num'] = int(proc_num/2)

  if ( 'lmp_openmp_num' in new_environ_dic ):
    lmp_openmp_num = new_environ_dic['lmp_openmp_num']
    if ( data_op.eval_str(lmp_openmp_num) == 1 ):
      new_environ_dic['lmp_openmp_num'] = int(lmp_openmp_num)
    else:
      log_info.log_error('Input error: lmp_openmp_num wrong, please check or set deepff/environ/lmp_openmp_num')
      exit()
  else:
    new_environ_dic['lmp_openmp_num'] = 2

  return new_deepmd_dic, new_lmp_dic, new_cp2k_dic, new_force_eval_dic, new_environ_dic
