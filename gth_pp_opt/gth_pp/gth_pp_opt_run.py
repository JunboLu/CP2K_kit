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
from CP2K_kit.gth_pp_opt.gth_pp import check_gth_opt
from CP2K_kit.gth_pp_opt.gth_pp import write_data
from CP2K_kit.gth_pp_opt.gth_pp import gen_atom_inp
from CP2K_kit.gth_pp_opt.gth_pp import init_step
from CP2K_kit.gth_pp_opt.gth_pp import step_reweight
from CP2K_kit.gth_pp_opt.gth_pp import weight_perturb
from CP2K_kit.gth_pp_opt.gth_pp import converg_perturb

#We know that this module is not beautiful enough, but we have to...

#Every generation begins as babies,
#Every baby must learn from scratch,
#Some generations reach high endpoints than others,
#Before they die off because of stupidity, dirtiness......

def get_min_step(process_dir):

  value = []
  step_index = []
  cmd = "ls | grep %s" % ('step')
  step_num = len(call.call_returns_shell(process_dir, cmd))
  for i in range(step_num):
    cmd = "grep %s step_%s/atom.out" %("'Final value of function'", str(i+1))
    return_tmp = call.call_returns_shell(process_dir, cmd)
    if ( return_tmp != [] ):
      return_tmp_split = data_op.split_str(return_tmp[0], ' ')
      if ( len(return_tmp_split) > 5 ):
        value.append(float(return_tmp_split[5]))
        step_index.append(i+1)

  min_value = min(value)
  min_value_index = value.index(min_value)
  min_step = step_index[min_value_index]

  return min_step

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
restart_index = gth_opt_param['restart_index']
micro_max_cycle = gth_opt_param['micro_max_cycle']
weight_1 = gth_opt_param['weight_1']
weight_2 = gth_opt_param['weight_2']
proc_1_func_conv = gth_opt_param['proc_1_func_conv']
proc_1_step_start = gth_opt_param['proc_1_step_start']

#Generate atom input file
element, val_elec_num, method = gen_atom_inp.gen_atom_inp(work_dir, gth_opt_param, weight_1)

#Get defined r_loc
line = linecache.getline(gth_pp_file, 3)
line_split = data_op.split_str(line, ' ')
r_loc_def = float(line_split[0])

#Run process_1
if ( restart_stage == 0 ):
  print ('Process_1: initial optimization', flush=True)
  init_step.gen_init_step(work_dir, gth_pp_file)

  start = proc_1_step_start
  end = start+parallel_num-1
  cycle = math.ceil((129-proc_1_step_start+1)/parallel_num)

  for i in range(cycle):
    init_step.run_init_step(work_dir, cp2k_exe, parallel_exe, element, val_elec_num, method, parallel_num, start, end)
    start = start+parallel_num
    end = end+parallel_num
    if ( end > 129 ):
      end = 129
  write_data.write_restart(work_dir, gth_opt_param, 1)

if ( restart_stage == 0 or restart_stage == 1 ):
  process_2_dir = ''.join((work_dir, '/process_2'))
  if ( os.path.exists(process_2_dir) ):
    restart_num = len(call.call_returns_shell(process_2_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_2_dir) and restart_num > 1 ):
    process_2_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_num-1)))
    min_step = get_min_step(process_2_restart_dir)
    gth_pp_file = ''.join((process_2_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  else:
    #Choose the lowest value of process_1
    step_index = []
    value = []
    r_loc = []
    wfn_state_1 = []
    process_1_dir = ''.join((work_dir, '/process_1'))
    for i in range(129):
      cmd = "grep %s step_%s/atom.out" %("'Final value of function'", str(i+1))
      return_func = call.call_returns_shell(process_1_dir, cmd)
      if ( len(return_func) != 0 and 'No such file' not in return_func[0] ):
        return_func_split = data_op.split_str(return_func[0], ' ')
        cmd = "grep %s step_%s/atom.out" %("'s-states N=    1'", str(i+1))
        return_wfn = call.call_returns_shell(process_1_dir, cmd)
        return_wfn_split = data_op.split_str(return_wfn[0], ' ')
        if ( len(return_func_split) > 5 ):
          value.append(float(return_func_split[5]))
          wfn_state_1.append(float(return_wfn_split[len(return_wfn_split)-2].strip('[')))
          step_index.append(i+1)
          opt_gth_pp_file = ''.join((process_1_dir, '/step_', str(i+1), '/GTH-PARAMETER'))
          line = linecache.getline(opt_gth_pp_file, 3)
          line_split = data_op.split_str(line, ' ')
          r_loc.append(float(line_split[0]))

    if ( len(step_index) != 129 ):
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
    for i in range(len(value)):
      if ( value[i] < proc_1_func_conv and abs(r_loc[i]-r_loc_def) < r_loc_conv ):
        value_proc.append(value[i])
        wfn_state_1_proc.append(wfn_state_1[i])
        step_index_proc.append(step_index[i])

    wfn_state_1_proc_abs = [abs(x) for x in wfn_state_1_proc]
    value_proc_asc, asc_order = data_op.get_list_order(value_proc, 'ascend', True)
    wfn_scale = [1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    for scale in wfn_scale:
      for i in asc_order:
        if ( wfn_state_1_proc_abs[i] <= scale*min(wfn_state_1_proc_abs) ):
          choosed_index = step_index_proc[i]
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

  #Run process_2
  print ('Process_2: automated step size optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_2_min_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_index)))
      min_step = get_min_step(process_2_min_restart_dir)
      gth_pp_file = ''.join((process_2_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = step_reweight.run_step_weight(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                  method, val_elec_num, python_exe, get_min_index, weight_1)
  write_data.write_restart(work_dir, gth_opt_param, 2, restart_index)

if ( restart_stage == 0 or restart_stage ==1 or restart_stage == 2 ):
  process_3_dir = ''.join((work_dir, '/process_3'))
  if ( os.path.exists(process_3_dir) ):
    restart_num = len(call.call_returns_shell(process_3_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_3_dir) and restart_num > 1 ):
    process_3_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_num-1)))
    min_step = get_min_step(process_3_restart_dir)
    gth_pp_file = ''.join((process_3_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  else:
    #Choose the lowest value of process_2
    process_2_min_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_index)))
    min_step = get_min_step(process_2_min_restart_dir)
    gth_pp_file = ''.join((process_2_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
    str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)

  #Run process_3
  print ('Process_3: weight pertubation optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_3_min_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_index)))
      min_step = get_min_step(process_3_min_restart_dir)
      gth_pp_file = ''.join((process_3_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = weight_perturb.run_weight_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                      method, val_elec_num, python_exe, get_min_index, weight_1)
  write_data.write_restart(work_dir, gth_opt_param, 3, restart_index)

if ( restart_stage == 0 or restart_stage ==1 or restart_stage == 2 or restart_stage == 3 ):
  process_4_dir = ''.join((work_dir, '/process_4'))
  if ( os.path.exists(process_4_dir) ):
    restart_num = len(call.call_returns_shell(process_4_dir, "ls | grep 'restart'"))
  else:
    restart_num = 0
  if ( os.path.exists(process_4_dir) and restart_num > 1 ):
    process_4_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_num-1)))
    min_step = get_min_step(process_4_restart_dir)
    gth_pp_file = ''.join((process_4_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  else:
    #Choose the lowest value of process_3
    process_3_min_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_index)))
    min_step = get_min_step(process_3_min_restart_dir)
    gth_pp_file = ''.join((process_3_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
    str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)

  #Run process_4
  print ('Process_4: convergence value pertubation optimization', flush=True)
  for i in range(micro_max_cycle):
    if ( i != 0 ):
      process_4_min_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_index)))
      min_step = get_min_step(process_4_min_restart_dir)
      gth_pp_file = ''.join((process_4_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))

    restart_index = converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, element, \
                                                        method, val_elec_num, python_exe, get_min_index, weight_2)

  process_4_min_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_index)))
  min_step = get_min_step(process_4_min_restart_dir)
  gth_pp_file = ''.join((process_4_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
  str_print = 'Success: choose the GTH paramter in %s' %(gth_pp_file)
  str_print = data_op.str_wrap(str_print, 80, '  ')
  print (str_print, flush=True)

  min_step_dir = ''.join((process_4_min_restart_dir, '/step_', str(min_step)))
  final_pp_dir = ''.join((work_dir, '/final'))
  call.call_simple_shell(work_dir, 'mkdir final')
  cmd = "cp %s %s" % (''.join((min_step_dir, '/*')), )
  call.call_simple_shell(min_step_dir, cmd)
