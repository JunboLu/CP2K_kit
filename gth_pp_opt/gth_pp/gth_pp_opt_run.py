#! /usr/env/bin python

import sys
import math
import copy
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import read_input
from CP2K_kit.tools import data_op
from CP2K_kit.gth_pp_opt.gth_pp import gen_atom_inp
from CP2K_kit.gth_pp_opt.gth_pp import init_step
from CP2K_kit.gth_pp_opt.gth_pp import step_reweight
from CP2K_kit.gth_pp_opt.gth_pp import weight_perturb
from CP2K_kit.gth_pp_opt.gth_pp import converg_perturb

#We know that this module is not beautiful enough, but we have to...

def get_min_step(process_dir):

  value = []
  step_index = []
  cmd = "ls | grep %s" % ('step')
  step_num = len(call.call_returns_shell(process_dir, cmd))
  for i in range(step_num):
    cmd = "grep %s step_%s/atom.out" %("'Final value of function'", str(i+1))
    return_temp = call.call_returns_shell(process_dir, cmd)
    if ( return_temp != [] ):
      return_temp_split = data_op.split_str(return_temp[0], ' ')
      if ( len(return_temp_split) > 5 ):
        value.append(return_temp_split[5])
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

job_type_param = read_input.dump_info(work_dir, inp_file, [run_type])

if ( 'cp2k_exe' in job_type_param[0].keys() ):
  cp2k_exe = job_type_param[0]['cp2k_exe']
else:
  print ('No cp2k executable file found, please set cp2k_exe')
  exit()

if ( 'parallel_exe' in job_type_param[0].keys() ):
  parallel_exe = job_type_param[0]['parallel_exe']
else:
  print ('No parallel executable file found, please set parallel_exe')
  exit()

if ( 'parallel_num' in job_type_param[0].keys() ):
  parallel_num = int(job_type_param[0]['parallel_num'])
else:
  parallel_num = 1

if ( 'init_gth_pp_file' in job_type_param[0].keys() ):
  gth_pp_file = job_type_param[0]['init_gth_pp_file']
else:
  print ('No file of initial guess gth parameter found, please set init_gth_pp_file')

#Generate atom input file
element, val_elec_num, method = gen_atom_inp.gen_atom_inp(work_dir, job_type_param[0])

#Run process_1
init_step.gen_init_step(work_dir, gth_pp_file)

start = 1
end = start+parallel_num-1
cycle = math.ceil(129/parallel_num)

for i in range(cycle):
  init_step.run_init_step(work_dir, cp2k_exe, parallel_exe, element, val_elec_num, method, parallel_num, start, end)
  start = start+parallel_num
  end = end+parallel_num
  if ( end > 129 ):
    end = 129

#Choose the lowest value of process_1
value = []
step_index = []
process_1_dir = ''.join((work_dir, '/process_1'))
for i in range(129):
  cmd = "grep %s step_%s/atom.out" %("'Final value of function'", str(i+1))
  return_temp = call.call_returns_shell(process_1_dir, cmd)
  if ( return_temp != [] ):
    return_temp_split = data_op.split_str(return_temp[0], ' ')
    if ( len(return_temp_split) > 5 ):
      value.append(return_temp_split[5])
      step_index.append(i+1)

if ( len(step_index) != 129 ):
  step_index_res = data_op.gen_list(1,129,1)
  step_index_res_copy = copy.deepcopy(step_index_res)
  for i in step_index:
    if i in step_index_res_copy:
      step_index_res.remove(i)
  print ('Warning: uncompleted steps are ', step_index_res)

min_value = min(value)
min_value_index = value.index(min_value)
min_step = step_index[min_value_index]

#Run process_2
gth_pp_file = ''.join((work_dir, '/process_1/step_', str(min_step), '/GTH-PARAMETER'))
restart_index = step_reweight.run_step_weight(work_dir, gth_pp_file, cp2k_exe, parallel_exe, \
                                              element, method, val_elec_num, python_exe, get_min_index)

#Choose the lowest value of process_2
process_2_min_restart_dir = ''.join((work_dir, '/process_2/restart', str(restart_index)))
min_step = get_min_step(process_2_min_restart_dir)

#Run process_3
gth_pp_file = ''.join((process_2_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
restart_index = weight_perturb.run_weight_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, \
                                                  element, method, val_elec_num, python_exe, get_min_index)

#Choose the lowest value of process_3
process_3_min_restart_dir = ''.join((work_dir, '/process_3/restart', str(restart_index)))
min_step = get_min_step(process_3_min_restart_dir)

#Run process_4
gth_pp_file = ''.join((process_3_min_restart_dir, '/step_', str(min_step), '/GTH-PARAMETER'))
restart_index = converg_perturb.run_converg_perturb(work_dir, gth_pp_file, cp2k_exe, parallel_exe, \
                                                    element, method, val_elec_num, python_exe, get_min_index)

process_4_min_restart_dir = ''.join((work_dir, '/process_4/restart', str(restart_index)))
min_step = get_min_step(process_4_min_restart_dir)

min_step_dir = ''.join((process_4_min_restart_dir, '/step_', str(min_step)))
final_pp_dir = ''.join((work_dir, '/final'))
call.call_simple_shell(work_dir, 'mkdir final')
cmd = "cp %s %s" % (''.join((min_step_dir, '/*')), )
call.call_simple_shell(min_step_dir, cmd)

