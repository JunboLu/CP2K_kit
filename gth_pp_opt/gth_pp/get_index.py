#! /usr/env/bin python

import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info

def get_index(file_dir):

  cmd = "ls | grep %s" % ('step_')
  step_num = len(call.call_returns_shell(file_dir, cmd))

  cmd = "grep %s %s/step_1/atom.inp" %("ELECTRON_CONFIGURATION", file_dir)
  elec_config_num = len(call.call_returns_shell(file_dir, cmd))

  step_index = []
  value = []
  vir_eigen_d = []

  for step in range(step_num):
    cmd = "grep %s %s/step_%s/atom.out" %("'Final value of function'", file_dir, str(step+1))
    return_func = call.call_returns_shell(file_dir, cmd)
    if ( len(return_func) != 0 and 'No such file' not in return_func[0] ):
      return_func_split = data_op.split_str(return_func[0], ' ')
      cmd = "grep %s %s/step_%s/atom.out" %('U1', file_dir, str(step+1))
      return_vir_u1_line = call.call_returns_shell(file_dir, cmd)
      return_vir_u1 = []
      for i in range(elec_config_num*4*2):
        return_vir_u1_split = data_op.split_str(return_vir_u1_line[i], ' ')
        return_vir_u1.append(float(return_vir_u1_split[5].split('[')[0]))
      vir_eigen_u1_d = 0.0
      for i in range(4*elec_config_num):
        if ( return_vir_u1[i+4*elec_config_num] < 1.0E-10 and return_vir_u1[i] < 1.0E-10 ):
          vir_eigen_u1_d = vir_eigen_u1_d - return_vir_u1[i+4*elec_config_num] + return_vir_u1[i]
        else:
          vir_eigen_u1_d = vir_eigen_u1_d + return_vir_u1[i+4*elec_config_num] - return_vir_u1[i]

      cmd = "grep %s %s/step_%s/atom.out" %('U2', file_dir, str(step+1))
      return_vir_u2_line = call.call_returns_shell(file_dir, cmd)
      return_vir_u2 = []
      for i in range(elec_config_num*4*2):
        return_vir_u2_split = data_op.split_str(return_vir_u2_line[i], ' ')
        return_vir_u2.append(float(return_vir_u2_split[5].split('[')[0]))
      vir_eigen_u2_d = 0.0
      for i in range(4*elec_config_num):
        if ( return_vir_u2[i+4*elec_config_num] < 1.0E-10 and return_vir_u2[i] < 1.0E-10 ):
          vir_eigen_u2_d = vir_eigen_u2_d - return_vir_u2[i+4*elec_config_num] + return_vir_u2[i]
        else:
          vir_eigen_u2_d = vir_eigen_u2_d + return_vir_u2[i+4*elec_config_num] - return_vir_u2[i]

      if ( len(return_func_split) > 5 ):
        value.append(float(return_func_split[5]))
        vir_eigen_d.append(vir_eigen_u1_d + vir_eigen_u2_d)
        step_index.append(step+1)

  vir_eigen_d_asc, asc_order = data_op.get_list_order(vir_eigen_d, 'ascend', True)
  scale_list = [1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18]
  for scale in scale_list:
    for i in asc_order:
      if ( value[i] <= scale*min(value) ):
        choosed_index = step_index[i]
        break
    if ( 'choosed_index' in locals() ):
      break

  if ( 'choosed_index' not in locals() ):
    log_info.log_error('Running error: no good parameter in %s' %(file_dir))
    exit()

  return choosed_index

if __name__ == '__main__':
  import sys
  file_name = str(sys.argv[1])
  m=get_index(file_name)
  print (m)
