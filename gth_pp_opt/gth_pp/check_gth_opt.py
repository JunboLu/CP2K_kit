#! /usr/env/bin python

import os
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op

def check_gth_opt(gth_opt_param):

  '''
  check_gth_pp: check the keywords of gth pp optimization.

  Args:
    gth_opt_param: dictionary
      gth_opt_param contains parameters for gth pp optimization.
  Returns:
    gth_opt_param: dictionary
      gth_opt_param is the revised gth_opt_param.
  '''

  valid_keyword = ['element', 'elec_config', 'elec_core_config', 'all_elec_method', 'xc_func', \
                   'relat_method', 'cp2k_exe', 'parallel_exe', 'init_gth_pp_file', 'restart_stage', \
                   'restart_index', 'micro_max_cycle', 'r_loc_conv', 'weight_1', 'weight_2', \
                   'proc_1_step_start', 'proc_1_func_conv', 'weight_pertub_1', 'weight_pertub_2', \
                   'weight_pertub_3', 'weight_pertub_4', 'consider_wfn_0']

  for key in gth_opt_param.keys():
    if key not in valid_keyword:
      log_info.log_error('Input error: %s is invalid keyword, please check the input file' %(key))
      exit()

  if ( 'element' not in gth_opt_param.keys() ):
    log_info.log_error('Input error: no element found, please set gth_pp_opt/element')
    exit()

  if ( 'elec_config' not in gth_opt_param.keys() ):
    log_info.log_error('Input error: no electron configuration found, please set gth_pp_opt/elec_config')
    exit()

  if ( 'elec_core_config' not in gth_opt_param.keys() ):
    log_info.log_error('Input error: no core electron configuration found, please set gth_pp_opt/elec_core_config')
    exit()

  if ( 'all_elec_method' in gth_opt_param.keys() ):
    all_elec_method = gth_opt_param['all_elec_method']
    if ( all_elec_method.upper() == 'HARTREE_FOCK' or all_elec_method.upper() == 'KOHN-SHAM' ):
      pass
    else:
      log_info.log_error('Input error: all_elec_method should be hartree-fock or kohn-sham, please check or reset gth_pp_opt/all_elec_method')
      exit()
  else:
    gth_opt_param['all_elec_method'] = 'HARTREE_FOCK'

  all_elec_method = gth_opt_param['all_elec_method']
  if ( all_elec_method.upper() == 'KOHN-SHAM' ):
    if ( 'xc_func' not in gth_opt_param.keys() ):
      gth_opt_param['xc_func'] = 'PBE'

  if ( 'relat_method' in gth_opt_param.keys() ):
    relat_method = gth_opt_param['relat_method']
    valid_relat_method = ['DKH(0)', 'DKH(1)', 'DKH(2)', 'DKH(3)', 'OFF', 'SCZORA(MP)', 'ZORA(MP)']
    if relat_method.upper() not in valid_relat_method:
      log_info.log_error('Input error: relativistic method is invalid, please check or reset gth_pp_opt/relat_method')
      exit()
  else:
    gth_opt_param['relat_method'] = 'DKH(3)'

  if ( 'consider_wfn_0' in gth_opt_param.keys() ):
    consider_wfn_0 = gth_opt_param['consider_wfn_0']
    if consider_wfn_0.upper() == 'TRUE' or consider_wfn_0.upper() == 'FALSE' :
      pass
    else:
      log_info.log_error('Input error: consider_wfn_0 should be bool, please check or reset gth_pp_opt/consider_wfn_0')
      exit()
  else:
    gth_opt_param['consider_wfn_0'] = 'true'

  if ( 'cp2k_exe' in gth_opt_param.keys() ):
    cp2k_exe = gth_opt_param['cp2k_exe']
    if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
      gth_opt_param['cp2k_exe'] = os.path.abspath(cp2k_exe)
    else:
      log_info.log_error('Input error: cp2k executable file does not exist, please check or reset gth_pp_opt/cp2k_exe')
      exit()
  else:
    log_info('Input error: no cp2k executable file found, please set gth_pp_opt/cp2k_exe')
    exit()

  if ( 'parallel_exe' in gth_opt_param.keys() ):
    parallel_exe = gth_opt_param['parallel_exe']
    if ( os.path.exists(os.path.abspath(parallel_exe)) ):
      gth_opt_param['parallel_exe'] = os.path.abspath(parallel_exe)
    else:
      log_info.log_error('Input error: parallel executable file does not exist, please check or reset gth_pp_opt/parallel_exe')
      exit()
  else:
    log_info('Input error: no parallel executable file found, please set gth_pp_opt/parallel_exe')
    exit()

  if ( 'init_gth_pp_file' in gth_opt_param.keys() ):
    init_gth_pp_file = gth_opt_param['init_gth_pp_file']
    if ( os.path.exists(os.path.abspath(init_gth_pp_file)) ):
      gth_opt_param['init_gth_pp_file'] = os.path.abspath(init_gth_pp_file)
    else:
      log_info.log_error('Input error: initial gth pp file does not exist, please check or reset gth_pp_opt/init_gth_pp_file')
      exit()
  else:
    log_info('Input error: no initial gth pp file found, please set gth_pp_opt/init_gth_pp_file')
    exit()

  if ( 'restart_stage' in gth_opt_param.keys() ):
    restart_stage = gth_opt_param['restart_stage']
    if ( data_op.eval_str(restart_stage) == 1 ):
      gth_opt_param['restart_stage'] = int(restart_stage)
    else:
      log_info.log_error('Input error: restart_stage should be integer, please check or reset gth_pp_opt/restart_stage')
      exit()
  else:
    gth_opt_param['restart_stage'] = 0

  if ( 'restart_index' in gth_opt_param.keys() ):
    restart_index = gth_opt_param['restart_index']
    if ( data_op.eval_str(restart_index) == 1 ):
      gth_opt_param['restart_index'] = int(restart_index)
    else:
      log_info.log_error('Input error: restart_index should be integer, please check or reset gth_pp_opt/restart_index')
      exit()
  else:
    gth_opt_param['restart_index'] = 0

  if ( 'proc_1_step_start' in gth_opt_param.keys() ):
    proc_1_step_start = gth_opt_param['proc_1_step_start']
    if ( data_op.eval_str(proc_1_step_start) == 1 ):
      gth_opt_param['proc_1_step_start'] = int(proc_1_step_start)
    else:
      log_info.log_error('Input error: proc_1_step_start should be integer, please check or reset gth_pp_opt/proc_1_step_start')
      exit()
  else:
    gth_opt_param['proc_1_step_start'] = 1

  if ( 'micro_max_cycle' in gth_opt_param.keys() ):
    micro_max_cycle = gth_opt_param['micro_max_cycle']
    if ( data_op.eval_str(micro_max_cycle) == 1 ):
      gth_opt_param['micro_max_cycle'] = int(micro_max_cycle)
    else:
      log_info.log_error('Input error: micro_max_cycle should be integer, please check or reset gth_pp_opt/micro_max_cycle')
      exit()
  else:
    gth_opt_param['micro_max_cycle'] = 1

  if ( 'r_loc_conv' in gth_opt_param.keys() ):
    r_loc_conv = gth_opt_param['r_loc_conv']
    if ( data_op.eval_str(r_loc_conv) == 1 or data_op.eval_str(r_loc_conv) == 2 ):
      gth_opt_param['r_loc_conv'] = float(r_loc_conv)
    else:
      log_info.log_error('Input error: r_loc_conv should be float or integer, please check or reset gth_pp_opt/r_loc_conv')
      exit()
  else:
    gth_opt_param['r_loc_conv'] = 0.005

  if ( 'proc_1_func_conv' in gth_opt_param.keys() ):
    proc_1_func_conv = gth_opt_param['proc_1_func_conv']
    if ( data_op.eval_str(proc_1_func_conv) == 1 or data_op.eval_str(proc_1_func_conv) == 2 ):
      gth_opt_param['proc_1_func_conv'] = float(proc_1_func_conv)
    else:
      log_info.log_error('Input error: proc_1_func_conv should be float or integer, please check or reset gth_pp_opt/proc_1_func_conv')
      exit()
  else:
    gth_opt_param['proc_1_func_conv'] = 10000.0

  if ( 'weight_1' in gth_opt_param.keys() ):
    weight_1 = gth_opt_param['weight_1']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_1):
      gth_opt_param['weight_1'] = [float(x) for x in weight_1]
    else:
      log_info.log_error('Input error: weight value in weight_1 should be integer or float, please check or reset gth_pp_opt/weight_1')
      exit()
  else:
    gth_opt_param['weight_1'] = [5.0, 2.0, 1.0]

  if ( 'weight_2' in gth_opt_param.keys() ):
    weight_2 = gth_opt_param['weight_2']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_2):
      gth_opt_param['weight_2'] = [float(x) for x in weight_2]
    else:
      log_info.log_error('Input error: weight value in weight_2 should be integer or float, please check or reset gth_pp_opt/weight_2')
      exit()
  else:
    gth_opt_param['weight_2'] = [2.0, 5.0, 1.0]

  if ( 'weight_pertub_1' in gth_opt_param.keys() ):
    weight_pertub_1 = gth_opt_param['weight_pertub_1']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_pertub_1):
      gth_opt_param['weight_pertub_1'] = [float(x) for x in weight_pertub_1]
    else:
      log_info.log_error('Input error: weight value in weight_pertub_1 should be integer or float, please check or reset gth_pp_opt/weight_pertub_1')
      exit()
  else:
    gth_opt_param['weight_pertub_1'] = [2.0, 5.0, 1.0]

  if ( 'weight_pertub_2' in gth_opt_param.keys() ):
    weight_pertub_2 = gth_opt_param['weight_pertub_2']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_pertub_2):
      gth_opt_param['weight_pertub_2'] = [float(x) for x in weight_pertub_2]
    else:
      log_info.log_error('Input error: weight value in weight_pertub_2 should be integer or float, please check or reset gth_pp_opt/weight_pertub_2')
      exit()
  else:
    gth_opt_param['weight_pertub_2'] = [5.0, 1.0, 2.0]

  if ( 'weight_pertub_3' in gth_opt_param.keys() ):
    weight_pertub_3 = gth_opt_param['weight_pertub_3']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_pertub_3):
      gth_opt_param['weight_pertub_3'] = [float(x) for x in weight_pertub_3]
    else:
      log_info.log_error('Input error: weight value in weight_pertub_3 should be integer or float, please check or reset gth_pp_opt/weight_pertub_3')
      exit()
  else:
    gth_opt_param['weight_pertub_3'] = [2.0, 1.0, 5.0]

  if ( 'weight_pertub_4' in gth_opt_param.keys() ):
    weight_pertub_4 = gth_opt_param['weight_pertub_4']
    if all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in weight_pertub_4):
      gth_opt_param['weight_pertub_4'] = [float(x) for x in weight_pertub_4]
    else:
      log_info.log_error('Input error: weight value in weight_pertub_4 should be integer or float, please check or reset gth_pp_opt/weight_pertub_4')
      exit()
  else:
    gth_opt_param['weight_pertub_4'] = [1.0, 5.0, 2.0]

  return gth_opt_param
