#! /usr/env/bin python

import os
import sys
import math
import copy
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import read_input
from CP2K_kit.tools import data_op
from CP2K_kit.deepff import sys_info
from CP2K_kit.gth_pp_opt.gth_pp import get_index
from CP2K_kit.gth_pp_opt.gth_pp import check_gth_opt
from CP2K_kit.gth_pp_opt.gth_pp import write_data
from CP2K_kit.gth_pp_opt.gth_pp import gen_atom_inp
from CP2K_kit.gth_pp_opt.gth_pp import init_step
from CP2K_kit.gth_pp_opt.gth_pp import step_reweight
from CP2K_kit.gth_pp_opt.gth_pp import weight_perturb
from CP2K_kit.gth_pp_opt.gth_pp import converg_perturb
from CP2K_kit.gth_pp_opt.gth_pp import elec_conf_perturb

#We know that this module is not beautiful enough, but we have to...

#Every generation begins as babies,
#Every baby must learn from scratch,
#Some generations reach high endpoints than others,
#Before they die off because of stupidity, dirtiness......

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
run_type = str(sys.argv[3])
python_exe = str(sys.argv[4])
get_min_index = ''.join((str(sys.argv[5]), '/gth_pp_opt/gth_pp/get_index.py'))

log_info.log_logo()

print (data_op.str_wrap('GTH_PP_OPT| PROGRAM STARTED IN %s' %(work_dir), 80), flush=True)
print ('GTH_PP_OPT| Input file name %s\n' %(inp_file), flush=True)
print ('GTH_PP_OPT'.center(80,'*'), flush=True)

parallel_num, proc_num_per_node, host, ssh = sys_info.get_host(work_dir)
if ( len(proc_num_per_node) != 1 ):
  log_info.log_error('The job is runing in %d nodes, it is better to use just one node' %(len(proc_num_per_node)), 'Warning')
  exit()

if ( parallel_num > 12 or parallel_num < 10 ):
  log_info.log_error('The job uses %d cores, it is better to use 10-12 cores' %(parallel_num), 'Warning')

job_type_param = read_input.dump_info(work_dir, inp_file, [run_type])
gth_opt_param = check_gth_opt.check_gth_opt(job_type_param[0])
cp2k_exe = gth_opt_param['cp2k_exe']
parallel_exe = gth_opt_param['parallel_exe']
gth_pp_file = gth_opt_param['init_gth_pp_file']
restart_stage = gth_opt_param['restart_stage']
r_loc_conv = gth_opt_param['r_loc_conv']
target_semi = gth_opt_param['target_semi']
target_val = gth_opt_param['target_val']
target_vir = gth_opt_param['target_vir']
micro_max_cycle = gth_opt_param['micro_max_cycle']
weight_1 = gth_opt_param['weight_1']
weight_2 = gth_opt_param['weight_2']
mix_max_cycle = gth_opt_param['mix_max_cycle']
mix_weight_1 = gth_opt_param['mix_weight_1']
mix_weight_2 = gth_opt_param['mix_weight_2']
weight_pertub_1 = gth_opt_param['weight_pertub_1']
weight_pertub_2 = gth_opt_param['weight_pertub_2']
weight_pertub_3 = gth_opt_param['weight_pertub_3']
weight_pertub_4 = gth_opt_param['weight_pertub_4']
converge_perturb_choice_1 = gth_opt_param['converge_perturb_choice_1']
converge_perturb_choice_2 = gth_opt_param['converge_perturb_choice_2']
converge_perturb_choice_3 = gth_opt_param['converge_perturb_choice_3']
converge_perturb_choice_4 = gth_opt_param['converge_perturb_choice_4']
converge_perturb_choice_5 = gth_opt_param['converge_perturb_choice_5']
converge_perturb_choice_6 = gth_opt_param['converge_perturb_choice_6']
converge_perturb_choice_7 = gth_opt_param['converge_perturb_choice_7']
converge_perturb_choice_8 = gth_opt_param['converge_perturb_choice_8']
converge_perturb_choice_9 = gth_opt_param['converge_perturb_choice_9']
converge_perturb_choice_10 = gth_opt_param['converge_perturb_choice_10']
converge_perturb_choice_11 = gth_opt_param['converge_perturb_choice_11']
converge_perturb_choice_12 = gth_opt_param['converge_perturb_choice_12']
elec_config_perturb_choice_1 = gth_opt_param['elec_config_perturb_choice_1']
elec_config_perturb_choice_2 = gth_opt_param['elec_config_perturb_choice_2']
elec_config_perturb_choice_3 = gth_opt_param['elec_config_perturb_choice_3']
elec_config_perturb_choice_4 = gth_opt_param['elec_config_perturb_choice_4']

proc_1_func_conv = gth_opt_param['proc_1_func_conv']
proc_1_step_start = gth_opt_param['proc_1_step_start']
weight_psir0 = gth_opt_param['weight_psir0']
weight_pot_node = gth_opt_param['weight_pot_node']
consider_wfn_0 = data_op.str_to_bool(gth_opt_param['consider_wfn_0'])
consider_r_loc = data_op.str_to_bool(gth_opt_param['consider_r_loc'])
opt_from_init = data_op.str_to_bool(gth_opt_param['opt_from_init'])

if ( 'elec_config_0' in gth_opt_param.keys()):
  elec_config = gth_opt_param['elec_config_0']
else:
  elec_config = gth_opt_param['elec_config']

#Generate atom input file
element, val_elec_num, method, elec_config_num = gen_atom_inp.gen_atom_inp(work_dir, gth_opt_param)

#Get defined r_loc
line = linecache.getline(gth_pp_file, 3)
line_split = data_op.split_str(line, ' ')
r_loc_def = float(line_split[0])

#Run process_1
if ( restart_stage == 0 ):
  print ('Process_1: initial optimization', flush=True)
  write_data.write_restart(work_dir, gth_opt_param, 0, gth_pp_file)
  init_step.gen_init_step(work_dir, gth_pp_file)

  start = proc_1_step_start
  for i in range(129):
    output_file = ''.join((work_dir, '/process_1/step_', str(i+1), '/atom.out'))
    cmd = "grep %s %s" %("'Final value of function'", output_file)
    return_func = call.call_returns_shell(''.join((work_dir, '/process_1')), cmd)
    gth_file = ''.join((work_dir, '/process_1/step_', str(i+1), '/GTH-PARAMETER'))
    if ( os.path.exists(gth_file) ):
      line_num = len(open(gth_file).readlines())
    else:
      line_num = 0
    if ( len(return_func) == 0 or 'No such file' in return_func[0] or line_num < 3 ):
      break
    else:
      start = i+1

  end = start+parallel_num-1
  if ( end > 129 ):
    end = 129
  cycle = math.ceil((129-start+1)/parallel_num)

  for i in range(cycle):
    init_step.run_init_step(work_dir, cp2k_exe, parallel_exe, element, val_elec_num, method, parallel_num, start, end)
    start = start+parallel_num
    end = end+parallel_num
    if ( end > 129 ):
      end = 129

  write_data.write_restart(work_dir, gth_opt_param, 1, gth_pp_file)

  #Choose the lowest value of process_1
if ( restart_stage == 0 or restart_stage == 1 ):
  step_index = []
  value = []
  r_loc = []
  wfn_state_1 = []
  eigen_dcharge_d = []
  process_1_dir = ''.join((work_dir, '/process_1'))
  for step in range(129):
    cmd = "grep %s %s/step_%s/atom.out" %("'Final value of function'", process_1_dir, str(step+1))
    return_func = call.call_returns_shell(process_1_dir, cmd)
    if ( len(return_func) != 0 and 'No such file' not in return_func[0] ):
      return_func_split = data_op.split_str(return_func[0], ' ')
      cmd = "grep %s %s/step_%s/atom.out" %("'s-states N=    1'", process_1_dir, str(step+1))
      return_wfn = call.call_returns_shell(process_1_dir, cmd)
      return_wfn_split = data_op.split_str(return_wfn[0], ' ')
      cmd = "grep %s %s/step_%s/atom.out" %('U1', process_1_dir, str(step+1))
      return_vir_u1_line = call.call_returns_shell(process_1_dir, cmd)
      return_vir_u1 = []
      for i in range(elec_config_num*4*2):
        return_vir_u1_split = data_op.split_str(return_vir_u1_line[i], ' ')
        return_vir_u1.append(float(return_vir_u1_split[5].split('[')[0]))
      vir_eigen_u1_d = 0.0
      for i in range(4*elec_config_num):
        if ( return_vir_u1[i+4*elec_config_num] < 1.0E-10 and return_vir_u1[i] < 1.0E-10 ):
          vir_eigen_u1_d = vir_eigen_u1_d - return_vir_u1[i+4*elec_config_num] + return_vir_u1[i]
        else:
          vir_eigen_u1_d = vir_eigen_u1_d + abs(return_vir_u1[i+4*elec_config_num]) - abs(return_vir_u1[i])

      cmd = "grep %s %s/step_%s/atom.out" %('U2', process_1_dir, str(step+1))
      return_vir_u2_line = call.call_returns_shell(process_1_dir, cmd)
      return_vir_u2 = []
      for i in range(elec_config_num*4*2):
        return_vir_u2_split = data_op.split_str(return_vir_u2_line[i], ' ')
        return_vir_u2.append(float(return_vir_u2_split[5].split('[')[0]))
      vir_eigen_u2_d = 0.0
      for i in range(4*elec_config_num):
        if ( return_vir_u2[i+4*elec_config_num] < 1.0E-10 and return_vir_u2[i] < 1.0E-10 ):
          vir_eigen_u2_d = vir_eigen_u2_d - return_vir_u2[i+4*elec_config_num] + return_vir_u2[i]
        else:
          vir_eigen_u2_d = vir_eigen_u2_d + abs(return_vir_u2[i+4*elec_config_num]) - abs(return_vir_u2[i])

      cmd = "grep %s %s/step_%s/atom.out" %("' SC '", process_1_dir, str(step+1))
      return_sc_line = call.call_returns_shell(process_1_dir, cmd)
      return_sc = []
      for i in range(len(return_sc_line)):
        return_sc_split_1 = data_op.split_str(return_sc_line[i], '[')
        return_sc_split = data_op.split_str(return_sc_split_1[1], ' ')
        return_sc.append(float(return_sc_split[1]))
      sc_dcharge_d = 0.0
      sc_num = int(len(return_sc)/(elec_config_num*2))
      for i in range(sc_num*elec_config_num):
        if ( return_sc[i+sc_num*elec_config_num] < 1.0E-10 and return_sc[i] < 1.0E-10 ):
          sc_dcharge_d = sc_dcharge_d - return_sc[i+sc_num*elec_config_num] + return_sc[i]
        else:
          sc_dcharge_d = sc_dcharge_d + abs(return_sc[i+sc_num*elec_config_num]) - abs(return_sc[i])

      cmd = "grep %s %s/step_%s/atom.out" %("' VA '", process_1_dir, str(step+1))
      return_val_line = call.call_returns_shell(process_1_dir, cmd)
      return_val = []
      for i in range(len(return_val_line)):
        return_val_split_1 = data_op.split_str(return_val_line[i], '[')
        return_val_split = data_op.split_str(return_val_split_1[1], ' ')
        return_val.append(float(return_val_split[1]))
      val_num = int(len(return_val)/(elec_config_num*2))
      val_dcharge_d = 0.0
      for i in range(val_num*elec_config_num):
        if ( return_val[i+val_num*elec_config_num] < 1.0E-10 and return_val[i] < 1.0E-10 ):
          val_dcharge_d = val_dcharge_d - return_val[i+val_num*elec_config_num] + return_val[i]
        else:
          val_dcharge_d = val_dcharge_d + abs(return_val[i+val_num*elec_config_num]) - abs(return_val[i])

      if ( len(return_func_split) > 5 ):
        value.append(float(return_func_split[5]))
        wfn_state_1.append(float(return_wfn_split[len(return_wfn_split)-2].strip('[')))
        eigen_dcharge_d.append(vir_eigen_u1_d+vir_eigen_u2_d+sc_dcharge_d+val_dcharge_d)
        step_index.append(step+1)
        opt_gth_pp_file = ''.join((process_1_dir, '/step_', str(step+1), '/GTH-PARAMETER'))
        line = linecache.getline(opt_gth_pp_file, 3)
        line_split = data_op.split_str(line, ' ')
        r_loc.append(float(line_split[0]))

  if ( len(step_index) == 0 ):
    log_info.log_error('Running error: no optimized paramter in %s, please check cp2k running in this directory' %(process_1_dir))
    exit()
  elif ( len(step_index) > 0 and len(step_index) < 129 ):
    step_index_res = data_op.gen_list(1,129,1)
    step_index_res_copy = copy.deepcopy(step_index_res)
    for i in step_index:
      if i in step_index_res_copy:
        step_index_res.remove(i)
    str_print = 'Warning: uncompleted steps are: %s' %(data_op.comb_list_2_str(step_index_res, ' '))
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)

  value_proc = []
  wfn_state_1_proc = []
  step_index_proc = []
  eigen_dcharge_d_proc = []
  if consider_r_loc:
    for i in range(len(value)):
      if ( value[i] < proc_1_func_conv and abs(r_loc[i]-r_loc_def) < r_loc_conv ):
        value_proc.append(value[i])
        wfn_state_1_proc.append(wfn_state_1[i])
        eigen_dcharge_d_proc.append(eigen_dcharge_d[i])
        step_index_proc.append(step_index[i])
  else:
    for i in range(len(value)):
      if ( value[i] < proc_1_func_conv ):
        value_proc.append(value[i])
        wfn_state_1_proc.append(wfn_state_1[i])
        eigen_dcharge_d_proc.append(eigen_d_charged[i])
        step_index_proc.append(step_index[i])

  if ( len(value_proc) == 0 ):
    log_info.log_error('Running error: no good parameters, maybe users need to reset proc_1_func_conv and r_loc_conv')
    exit()

  eigen_dcharge_d_asc, asc_order = data_op.get_list_order(eigen_dcharge_d_proc, 'ascend', True)
  value_proc_asc = data_op.get_list_order(value_proc, 'ascend')
  value_scale_init = float(value_proc_asc[1]/value_proc_asc[0])
  wfn_state_1_proc_abs = [abs(x) for x in wfn_state_1_proc]
  wfn_scale_list = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  value_scale_list = []
  for i in range(21):
    value_scale_list.append(value_scale_init+0.02*(i+1))

  for i in range(21):
    for j in asc_order:
      if consider_wfn_0:
        if ( wfn_state_1_proc_abs[j] <= wfn_scale_list[i]*min(wfn_state_1_proc_abs) and \
             value_proc[j] <= value_scale_list[i]*min(value_proc) ):
          choosed_index = step_index_proc[j]
          break
      else:
        if ( value_proc[j] <= value_scale_list[i]*min(value_proc) ):
          choosed_index = step_index_proc[j]
          break
    if ( 'choosed_index' in locals() ):
      break

  if ( 'choosed_index' not in locals() ):
    log_info.log_error('Running error: no good parameter in process_1')
    exit()
  gth_pp_file = ''.join((work_dir, '/process_1/step_', str(choosed_index), '/GTH-PARAMETER'))
  str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
  str_print = data_op.str_wrap(str_print, 80, '  ')
  print (str_print, flush=True)

  write_data.write_restart(work_dir, gth_opt_param, 2, gth_pp_file)

if ( restart_stage == 0 or restart_stage == 1 or restart_stage == 2 ):
  process_2_dir = ''.join((work_dir, '/process_2'))
  if ( os.path.exists(process_2_dir) ):
    restart_num = len(call.call_returns_shell(process_2_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_2_dir) and restart_num > 1 and not opt_from_init ):
    process_2_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_num-1)))
    min_step = get_index(process_2_restart_dir)
    gth_pp_file = ''.join((process_2_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

  #Run process_2
  print ('Process_2: automated step size optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_2_min_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_index)))
      min_step = get_index(process_2_min_restart_dir)
      gth_pp_file = ''.join((process_2_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = step_reweight.run_step_weight(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                  method, val_elec_num, python_exe, get_min_index, weight_1, \
                                                  target_semi, target_val, target_vir)

  #Choose the lowest value of process_2
  process_2_min_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_index)))
  min_step = get_index(process_2_min_restart_dir)
  gth_pp_file = ''.join((process_2_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
  str_print = data_op.str_wrap(str_print, 80, '  ')
  print (str_print, flush=True)

  write_data.write_restart(work_dir, gth_opt_param, 3, gth_pp_file)

if ( restart_stage == 0 or restart_stage ==1 or restart_stage == 2 or restart_stage == 3 ):
  process_3_dir = ''.join((work_dir, '/process_3'))
  if ( os.path.exists(process_3_dir) ):
    restart_num = len(call.call_returns_shell(process_3_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_3_dir) and restart_num > 1 and not opt_from_init ):
    process_3_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_num-1)))
    min_step = get_index(process_3_restart_dir)
    gth_pp_file = ''.join((process_3_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

  #Run process_3
  print ('Process_3: weight pertubation optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_3_min_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_index)))
      min_step = get_index(process_3_min_restart_dir)
      gth_pp_file = ''.join((process_3_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = weight_perturb.run_weight_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                      method, val_elec_num, python_exe, get_min_index, weight_1, \
                                                      weight_pertub_1, weight_pertub_2, weight_pertub_3, weight_pertub_4, \
                                                      target_semi, target_val, target_vir)

  #Choose the lowest value of process_3
  process_3_min_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_index)))
  min_step = get_index(process_3_min_restart_dir)
  gth_pp_file = ''.join((process_3_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
  str_print = data_op.str_wrap(str_print, 80, '  ')
  print (str_print, flush=True)

  write_data.write_restart(work_dir, gth_opt_param, 4, gth_pp_file)

if ( restart_stage == 0 or restart_stage ==1 or restart_stage == 2 or restart_stage == 3 or restart_stage == 4 ):
  process_4_dir = ''.join((work_dir, '/process_4'))
  if ( os.path.exists(process_4_dir) ):
    restart_num = len(call.call_returns_shell(process_4_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_4_dir) and restart_num > 1 and not opt_from_init ):
    process_4_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_num-1)))
    min_step = get_index(process_4_restart_dir)
    gth_pp_file = ''.join((process_4_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

  #Run process_4
  print ('Process_4: convergence value pertubation optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_4_min_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_index)))
      min_step = get_index(process_4_min_restart_dir)
      gth_pp_file = ''.join((process_4_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                        method, val_elec_num, python_exe, get_min_index, weight_2, \
                                                        converge_perturb_choice_1, converge_perturb_choice_2, \
                                                        converge_perturb_choice_3, converge_perturb_choice_4, \
                                                        converge_perturb_choice_5, converge_perturb_choice_6, \
                                                        converge_perturb_choice_7, converge_perturb_choice_8, \
                                                        converge_perturb_choice_9, converge_perturb_choice_10, \
                                                        converge_perturb_choice_11, converge_perturb_choice_12, \
                                                        target_semi, target_val, target_vir)

if ( restart_stage == 5 ):

  elec_config_perturb = []
  elec_config_str = data_op.comb_list_2_str(elec_config, ' ')
  elec_config_perturb_choice_1_str = data_op.comb_list_2_str(elec_config_perturb_choice_1, ' ')
  elec_config_perturb_choice_2_str = data_op.comb_list_2_str(elec_config_perturb_choice_2, ' ')
  elec_config_perturb_choice_3_str = data_op.comb_list_2_str(elec_config_perturb_choice_3, ' ')
  elec_config_perturb_choice_4_str = data_op.comb_list_2_str(elec_config_perturb_choice_4, ' ')
  elec_config_perturb.append(elec_config_perturb_choice_1)
  elec_config_perturb.append(elec_config_perturb_choice_2)
  elec_config_perturb.append(elec_config_perturb_choice_3)
  elec_config_perturb.append(elec_config_perturb_choice_4)
  choice_num = 0
  for i in range(4):
    if ( elec_config_perturb[i] != 'none' ):
      choice_num = choice_num+1

  if ( choice_num >= 1 ):
    restart_index = elec_conf_perturb.run_elec_config_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                                method, val_elec_num, python_exe, get_min_index, choice_num, elec_config_str, \
                                                                elec_config_perturb_choice_1_str, elec_config_perturb_choice_2_str, \
                                                                elec_config_perturb_choice_3_str, elec_config_perturb_choice_4_str, \
                                                                target_semi, target_val, target_vir)


if ( restart_stage == 6 ):
  for i in range(100):
    #Run process_4
    restart_index = converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                        method, val_elec_num, python_exe, get_min_index, mix_weight_1, \
                                                        converge_perturb_choice_1, converge_perturb_choice_2, \
                                                        converge_perturb_choice_3, converge_perturb_choice_4, \
                                                        converge_perturb_choice_5, converge_perturb_choice_6, \
                                                        converge_perturb_choice_7, converge_perturb_choice_8, \
                                                        converge_perturb_choice_9, converge_perturb_choice_10, \
                                                        converge_perturb_choice_11, converge_perturb_choice_12, \
                                                        target_semi, target_val, target_vir, mix_max_cycle)
    min_restart_dir = ''.join((work_dir, '/process_mix/restart', str(restart_index)))
    min_step = get_index(min_restart_dir)
    gth_pp_file = ''.join((min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
    #Run process_3
    restart_index = weight_perturb.run_weight_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                      method, val_elec_num, python_exe, get_min_index, mix_weight_2, \
                                                      mix_weight_1, weight_pertub_2, weight_pertub_3, weight_pertub_4, \
                                                      target_semi, target_val, target_vir, mix_max_cycle)
    min_restart_dir = ''.join((work_dir, '/process_mix/restart', str(restart_index)))
    min_step = get_index(min_restart_dir)
    gth_pp_file = ''.join((min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
    #Run process_4
    restart_index = converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                        method, val_elec_num, python_exe, get_min_index, mix_weight_2,
                                                        converge_perturb_choice_1, converge_perturb_choice_2, \
                                                        converge_perturb_choice_3, converge_perturb_choice_4, \
                                                        converge_perturb_choice_5, converge_perturb_choice_6, \
                                                        converge_perturb_choice_7, converge_perturb_choice_8, \
                                                        converge_perturb_choice_9, converge_perturb_choice_10, \
                                                        converge_perturb_choice_11, converge_perturb_choice_12, \
                                                        target_semi, target_val, target_vir, mix_max_cycle)

    min_restart_dir = ''.join((work_dir, '/process_mix/restart', str(restart_index)))
    min_step = get_index(min_restart_dir)
    gth_pp_file = ''.join((min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
    #Run process_3
    restart_index = weight_perturb.run_weight_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                      method, val_elec_num, python_exe, get_min_index, mix_weight_1, \
                                                      mix_weight_2, weight_pertub_2, weight_pertub_3, weight_pertub_4, \
                                                      target_semi, target_val, target_vir, mix_max_cycle)

    min_restart_dir = ''.join((work_dir, '/process_mix/restart', str(restart_index)))
    min_step = get_index(min_restart_dir)
    gth_pp_file = ''.join((min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

  process_4_min_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_index)))
  min_step = get_index(process_4_min_restart_dir)
  gth_pp_file = ''.join((process_4_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
  str_print = data_op.str_wrap(str_print, 80, '  ')
  print (str_print, flush=True)

  min_step_dir = ''.join((process_4_min_restart_dir, '/step_', str(min_step)))
  final_pp_dir = ''.join((work_dir, '/final'))
  call.call_simple_shell(work_dir, 'mkdir final')
  cmd = "cp %s %s" % (''.join((min_step_dir, '/*')), final_pp_dir)
  call.call_simple_shell(min_step_dir, cmd)
