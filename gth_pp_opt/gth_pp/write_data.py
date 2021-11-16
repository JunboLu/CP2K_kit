#! /usr/env/bin python

from CP2K_kit.tools import data_op

def write_restart(work_dir, gth_pp_opt_dic, restart_stage):

  '''
  write_restart: write restarting input file.

  Args:
    work_dir: string
      work_dir is the working directory.
    gth_pp_opt_dic: dictionary
      gth_pp_opt_dic contains keywords used to generate atom input file.
    restart_stage: int
      restart_stage is the stage of restarting
  Returns:
    none
  '''

  restart_file_name = ''.join((work_dir, '/CP2K_KIT.restart'))
  restart_file = open(restart_file_name, 'w')

  restart_file.write('&global\n')
  restart_file.write('  run_type gth_pp_opt\n')
  restart_file.write('&end global\n')
  restart_file.write('\n')
  restart_file.write('&gth_pp_opt\n')

  for key in gth_pp_opt_dic.keys():
    if ( key != 'restart_stage' and key != 'elec_config' and key != 'elec_core_config' ):
      restart_file.write('  %s %s\n' %(key, gth_pp_opt_dic[key]))
    if ( key == 'elec_config' or key == 'elec_core_config' ):
      key_str = data_op.comb_list_2_str(gth_pp_opt_dic[key], ' ')
      restart_file.write('  %s %s\n' %(key, key_str))
  restart_file.write('  restart_stage %d\n' %(restart_stage))
  restart_file.write('&end gth_pp_opt\n')
  restart_file.close()