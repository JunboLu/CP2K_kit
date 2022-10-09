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
                   'max_interval', 'md_steps', 'restart_step', 'colvar_id']

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
    if ( os.path.exists(os.path.abspath(os.path.expanduser(cp2k_inp_file))) ):
      cmd_dic['cp2k_inp_file'] = os.path.abspath(os.path.expanduser(cp2k_inp_file))
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

  if ( 'max_interval' in cmd_dic.keys() ):
    max_interval = cmd_dic['max_interval']
    if ( data_op.eval_str(max_interval) == 1 ):
      cmd_dic['max_interval'] = int(max_interval)
    else:
      log_info.log_error('Input error: macro steps should be integer, please check or reset thermo_int/constraint_md/max_interval')
      exit()
  else:
    log_info.log_error('Input error: no macro steps, please set thermo_int/constraint_md/max_interval')
    exit()

  if ( 'md_steps' in cmd_dic.keys() ):
    md_steps = cmd_dic['md_steps']
    if ( data_op.eval_str(md_steps) == 1 ):
      cmd_dic['md_steps'] = int(md_steps)
    else:
      log_info.log_error('Input error: md steps should be integer, please check or reset thermo_int/constraint_md/md_steps')
      exit()
  else:
    log_info.log_error('Input error: no md steps, please set thermo_int/constraint_md/md_steps')
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

  mix_dic_valid_key = ['mix_inp_file', 'max_interval', 'md_steps', 'restart_step']

  for key in mix_dic.keys():
    if key not in mix_dic_valid_key:
      log_info.log_error('Input error: %s is invalid keyword, please check or reset thermo_int/mix_force_eval')
      exit()

  if ( 'mix_inp_file' in mix_dic.keys() ):
    mix_inp_file = mix_dic['mix_inp_file']
    if ( os.path.exists(os.path.abspath(os.path.expanduser(mix_inp_file))) ):
      mix_dic['mix_inp_file'] = os.path.abspath(os.path.expanduser(mix_inp_file))
    else:
      log_info.log_error('Input error: %s does not exist' %(mix_inp_file))
      exit()
  else:
    log_info.log_error('Input error: no cp2k mixing force_eval input file, please set thermo_int/mix_force_eval/mix_inp_file')
    exit()

  if ( 'max_interval' in mix_dic.keys() ):
    max_interval = mix_dic['max_interval']
    if ( data_op.eval_str(max_interval) == 1 ):
      mix_dic['max_interval'] = int(max_interval)
    else:
      log_info.log_error('Input error: macro steps should be integer, please check or reset thermo_int/mix_force_eval/max_interval')
      exit()
  else:
    mix_dic['max_interval'] = 10000

  if ( 'md_steps' in mix_dic.keys() ):
    md_steps = mix_dic['md_steps']
    if ( data_op.eval_str(md_steps) == 1 ):
      mix_dic['md_steps'] = int(md_steps)
    else:
      log_info.log_error('Input error: md steps should be integer, please check or reset thermo_int/mix_force_eval/md_steps')
      exit()
  else:
    mix_dic['md_steps'] = 1

  if ( 'restart_step' in mix_dic.keys() ):
    restart_step = mix_dic['restart_step']
    if ( data_op.eval_str(restart_step) == 1 ):
      mix_dic['restart_step'] = int(restart_step)
    else:
      log_info.log_error('Input error: restart step should be integer, please check or reset thermo_int/mix_force_eval/restart_step')
      exit()
  else:
    mix_dic['restart_step'] = 1

  return mix_dic
