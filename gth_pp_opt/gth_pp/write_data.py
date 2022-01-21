#! /usr/env/bin python

from CP2K_kit.tools import data_op

def write_restart(work_dir, gth_pp_opt_dic, restart_stage, gth_pp_file):

  '''
  write_restart: write restarting input file.

  Args:
    work_dir: string
      work_dir is the working directory.
    gth_pp_opt_dic: dictionary
      gth_pp_opt_dic contains keywords used to generate atom input file.
    restart_stage: int
      restart_stage is the stage of restarting
    gth_pp_file: string
      gth_pp_file is the name of gth pp file.
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
    if ( key != 'restart_stage' and \
         key != 'elec_config' and \
         key != 'elec_core_config' and \
         key != 'restart_index' and \
         key != 'weight_1' and \
         key != 'weight_2' and \
         key != 'weight_pertub_1' and \
         key != 'weight_pertub_2' and \
         key != 'weight_pertub_3' and \
         key != 'weight_pertub_4' and \
         key != 'init_gth_pp_file' ):
      restart_file.write('  %s %s\n' %(key, gth_pp_opt_dic[key]))
    if ( key == 'elec_config' or key == 'elec_core_config' ):
      if ( isinstance(gth_pp_opt_dic[key], list) ):
        key_str = data_op.comb_list_2_str(gth_pp_opt_dic[key], ' ')
        restart_file.write('  %s %s\n' %(key, key_str))
      else:
        restart_file.write('  %s %s\n' %(key, gth_pp_opt_dic[key]))
    if ( key == 'weight_1' or \
         key == 'weight_2' or \
         key == 'weight_pertub_1' or \
         key == 'weight_pertub_2' or \
         key == 'weight_pertub_3' or \
         key == 'weight_pertub_4'):
      key_str = data_op.comb_list_2_str(gth_pp_opt_dic[key], ' ')
      restart_file.write('  %s %s\n' %(key, key_str))
  restart_file.write('  init_gth_pp_file %s\n' %(gth_pp_file))
  restart_file.write('  restart_stage %d\n' %(restart_stage))
  restart_file.write('&end gth_pp_opt\n')
  restart_file.close()
