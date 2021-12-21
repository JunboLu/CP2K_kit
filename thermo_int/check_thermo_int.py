#! /usr/env/bin python

import os
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op

def check_constraint_md_inp(cmd_dic):

  '''
  check_constraint_md_inp: check the input of constraint md.

  Args:
    cmd_dic: dictionary
      cmd_dic contains the parameter for constraint md.
  Returns:
    cmd_dic: dictionary
      cmd_dic is the revised cmd_dic.
  '''

  cmd_valid_key = ['run_type', 'cp2k_inp_file', 'init_value', 'end_value', \
                   'macro_steps', 'micro_steps', 'restart_step', 'colvar_id', \
                   'cp2k_exe', 'cp2k_mpi_num']
  for key in cmd_dic.keys():
    if key not in cmd_valid_key:
      log_info.log_error('Input error: %s is invalid keyword, please check or reset thermo_int/constraint_md')
      exit()

  if ( 'run_type' in cmd_dic.keys() ):
    run_type = cmd_dic['run_type']
    if ( run_type == 'slow_growth' or run_type == 'finite_points' ):
      pass
    else:
      log_info.log_error('Input error: only slow_growth and finite_points are supported for run_type, please check or reset thermo_int/constraint_md/run_type')
      exit()
  else:
    log_info.log_error('Input error: no run_type, please set thermo_int/constraint_md/run_type')
    exit()

  if ( 'cp2k_inp_file' in cmd_dic.keys() ):
    cp2k_inp_file = cmd_dic['cp2k_inp_file']
    if ( os.path.exists(os.path.abspath(cp2k_inp_file)) ):
      cmd_dic['cp2k_inp_file'] = os.path.abspath(cp2k_inp_file)
    else:
      log_info.log_error('Input error: %s does not exist' %(cp2k_inp_file))
      exit()
  else:
    log_info.log_error('Input error: no cp2k input file, please set thermo_int/constraint_md/cp2k_inp_file')
    exit()

  if ( 'init_value' in cmd_dic.keys() ):
    init_value = cmd_dic['init_value']
    if ( data_op.eval_str(init_value) == 1 or data_op.eval_str(init_value) == 2 ):
      cmd_dic['init_value'] = float(init_value)
    else:
      log_info.log_error('Input error: initial value should be integer or float, please check or reset thermo_int/constraint_md/init_value')
      exit()
  else:
    log_info.log_error('Input error: no initial value, please set thermo_int/constraint_md/init_value')
    exit()

  if ( 'end_value' in cmd_dic.keys() ):
    end_value = cmd_dic['end_value']
    if ( data_op.eval_str(end_value) == 1 or data_op.eval_str(end_value) == 2 ):
      cmd_dic['end_value'] = float(end_value)
    else:
      log_info.log_error('Input error: endding value should be integer or float, please check or reset thermo_int/constraint_md/end_value')
      exit()
  else:
    log_info.log_error('Input error: no endding value, please set thermo_int/constraint_md/end_value')
    exit()

  if ( 'macro_steps' in cmd_dic.keys() ):
    macro_steps = cmd_dic['macro_steps']
    if ( data_op.eval_str(macro_steps) == 1 ):
      cmd_dic['macro_steps'] = int(macro_steps)
    else:
      log_info.log_error('Input error: macro steps should be integer, please check or reset thermo_int/constraint_md/macro_steps')
      exit()
  else:
    log_info.log_error('Input error: no macro steps, please set thermo_int/constraint_md/macro_steps')
    exit()

  if ( 'micro_steps' in cmd_dic.keys() ):
    micro_steps = cmd_dic['micro_steps']
    if ( data_op.eval_str(micro_steps) == 1 ):
      cmd_dic['micro_steps'] = int(micro_steps)
    else:
      log_info.log_error('Input error: micro steps should be integer, please check or reset thermo_int/constraint_md/micro_steps')
      exit()
  else:
    log_info.log_error('Input error: no micro steps, please set thermo_int/constraint_md/micro_steps')
    exit()

  if ( 'restart_step' in cmd_dic.keys() ):
    restart_step = cmd_dic['restart_step']
    if ( data_op.eval_str(restart_step) == 1 ):
      cmd_dic['restart_step'] = int(restart_step)
    else:
      log_info.log_error('Input error: restart step should be integer, please check or reset thermo_int/constraint_md/restart_step')
      exit()
  else:
    cmd_dic['restart_step'] = 1

  if ( 'colvar_id' in cmd_dic.keys() ):
    colvar_id = cmd_dic['colvar_id']
    if ( data_op.eval_str(colvar_id) == 1 ):
      cmd_dic['colvar_id'] = int(colvar_id)
    else:
      log_info.log_error('Input error: colvar id should be integer, please check or reset thermo_int/constraint_md/colvar_id')
      exit()
  else:
    cmd_dic['colvar_id'] = 1

  if ( 'cp2k_exe' in cmd_dic.keys() ):
    cp2k_exe = cmd_dic['cp2k_exe']
    if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
      cmd_dic['cp2k_exe'] = os.path.abspath(cp2k_exe)
    else:
      log_info.log_error('Input error: %s does not exists')
      exit()
  else:
    log_info.log_error('Input error: no cp2k executable file, please set thermo_int/constraint_md/cp2k_exe')
    exit()

  if ( 'cp2k_mpi_num' in cmd_dic.keys() ):
    cp2k_mpi_num = cmd_dic['cp2k_mpi_num']
    if ( data_op.eval_str(cp2k_mpi_num) == 1 ):
      cmd_dic['cp2k_mpi_num'] = int(cp2k_mpi_num)
    else:
      log_info.log_error('Input error: cp2k mpi number should be integer, please check or reset thermo_int/constraint_md/cp2k_mpi_num')
      exit()
  else:
    log_info.log_error('Input error: no cp2k mpi number, please set thermo_int/constraint_md/cp2k_mpi_num')
    exit()

  return cmd_dic

def check_mix_force_eval_inp(mix_dic):

  '''
  check_mix_force_eval_inp: check the input of mixing force eval.

  Args:
    mix_dic: dictionary
      mix_dic contains the parameter for mixing force eval.
  Returns:
    mix_dic: dictionary
      mix_dic is the revised mix_dic.
  '''

  mix_dic_valid_key = ['mix_inp_file', 'macro_steps', 'micro_steps', 'restart_step', 'cp2k_exe']
  for key in mix_dic.keys():
    if key not in mix_dic_valid_key:
      log_info.log_error('Input error: %s is invalid keyword, please check or reset thermo_int/mix_force_eval')
      exit()

  if ( 'mix_inp_file' in mix_dic.keys() ):
    mix_inp_file = mix_dic['mix_inp_file']
    if ( os.path.exists(os.path.abspath(mix_inp_file)) ):
      mix_dic['mix_inp_file'] = os.path.abspath(mix_inp_file)
    else:
      log_info.log_error('Input error: %s does not exist' %(mix_inp_file))
      exit()
  else:
    log_info.log_error('Input error: no cp2k mixing force_eval input file, please set thermo_int/mix_force_eval/mix_inp_file')
    exit()

  if ( 'macro_steps' in mix_dic.keys() ):
    macro_steps = mix_dic['macro_steps']
    if ( data_op.eval_str(macro_steps) == 1 ):
      mix_dic['macro_steps'] = int(macro_steps)
    else:
      log_info.log_error('Input error: macro steps should be integer, please check or reset thermo_int/mix_force_eval/macro_steps')
      exit()
  else:
    mix_dic['macro_steps'] = 10000

  if ( 'micro_steps' in mix_dic.keys() ):
    micro_steps = mix_dic['micro_steps']
    if ( data_op.eval_str(micro_steps) == 1 ):
      mix_dic['micro_steps'] = int(micro_steps)
    else:
      log_info.log_error('Input error: micro steps should be integer, please check or reset thermo_int/mix_force_eval/micro_steps')
      exit()
  else:
    mix_dic['micro_steps'] = 1

  if ( 'restart_step' in mix_dic.keys() ):
    restart_step = mix_dic['restart_step']
    if ( data_op.eval_str(restart_step) == 1 ):
      mix_dic['restart_step'] = int(restart_step)
    else:
      log_info.log_error('Input error: restart step should be integer, please check or reset thermo_int/mix_force_eval/restart_step')
      exit()
  else:
    mix_dic['restart_step'] = 1

  if ( 'cp2k_exe' in mix_dic.keys() ):
    cp2k_exe = mix_dic['cp2k_exe']
    if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
      mix_dic['cp2k_exe'] = os.path.abspath(cp2k_exe)
    else:
      log_info.log_error('Input error: %s does not exists')
      exit()
  else:
    log_info.log_error('Input error: no cp2k executable file, please set thermo_int/mix_force_eval/cp2k_exe')
    exit()

  return mix_dic
