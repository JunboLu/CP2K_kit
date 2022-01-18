#! /usr/env/bin python

import os
import linecache
import subprocess
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import file_tools
from CP2K_kit.thermo_int import check_thermo_int

def mix_force_eval(mix_inp_file, max_interval, md_steps, restart_step, cp2k_exe, work_dir):

  '''
  mix_force_eval: mixing force eval calculation

  Args:
    mix_inp_file: string
      mix_inp_file is the cp2k mixing force eval file.
    max_interval: int
      max_interval is the totals steps.
    md_steps: int
      md_steps is the steps for each md.
    restart_step: int
      restart_step is the restaring step.
    cp2k_exe: string
      cp2k_exe is the cp2k executable file.
    work_dir: string
      work_dir is the workding directory.
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(mix_inp_file, work_dir)

  #Get project name
  line_num = file_tools.grep_line_num("'PROJECT'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    proj_line_num = line_num[0]
  else:
    line_num = file_tools.grep_line_num("'PROJECT_NAME'", upper_file_name_abs, work_dir)
    if ( line_num != 0 ):
      proj_line_num = line_num[0]
    else:
      log_info.log_error('Input error: no project name in %s file.' %(mix_inp_file))
      exit()

  cmd = "rm %s" %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

  proj_line = linecache.getline(mix_inp_file, proj_line_num)
  proj_line_split = data_op.split_str(proj_line, ' ', '\n')
  proj_name = proj_line_split[len(proj_line_split)-1]

  restart_file_name = ''.join((proj_name, '-1.restart'))
  restart_file_name_abs = ''.join((work_dir, '/', restart_file_name))
  mix_ene_file_name = ''.join((proj_name, '-mix-1.ener'))
  mix_ene_file_name_abs = ''.join((work_dir, '/', mix_ene_file_name))

  if ( restart_step == 1 ):
    cmd = "cp %s %s" %(mix_inp_file, restart_file_name)
    call.call_simple_shell(work_dir, cmd)

  incre_value = float(1.0/max_interval)

  m = 0
  while True:
    inp_trans_file_name = ''.join(('input_mix', str(m), '.inp'))
    inp_trans_file_name_abs = ''.join((work_dir, '/', inp_trans_file_name))
    if ( os.path.exists(inp_trans_file_name_abs) ):
      m = m+1
    else:
      break

  for i in range(restart_step, max_interval+2, 1):
    #Write cp2k_kit restart file
    while True:
      restart_inp_file_name_abs = ''.join((work_dir, '/CP2K_KIT.restart'))
      if ( os.path.exists(restart_inp_file_name_abs) ):
        cmd = 'rm %s' %('CP2K_KIT.restart')
        call.call_simple_shell(work_dir, cmd)
      else:
        break

    if ( i == restart_step and restart_step != 1 ):
      cmd = "cp %s %s" %('RESTART_BAK_FILE', restart_file_name)
      call.call_simple_shell(work_dir, cmd)

    if ( os.path.exists(mix_ene_file_name_abs) ):
      whole_line_num = len(open(mix_ene_file_name_abs, 'r').readlines())
      if ( whole_line_num > 2*(restart_step-1)*md_steps ):
        cmd = "sed -i '%d,%dd' %s" %((restart_step-1)*md_steps+1, whole_line_num, mix_ene_file_name)
        call.call_simple_shell(work_dir, cmd)

    restart_inp_file = open(restart_inp_file_name_abs, 'w')
    restart_inp_file.write('&global\n')
    restart_inp_file.write('  run_type thermo_int\n')
    restart_inp_file.write('  thermo_type mix_force_eval\n')
    restart_inp_file.write('&end global\n')
    restart_inp_file.write('\n')

    restart_inp_file.write('&thermo_int\n')
    restart_inp_file.write('  &mix_force_eval\n')
    restart_inp_file.write('    mix_inp_file %s\n' %(mix_inp_file))
    restart_inp_file.write('    max_interval %d\n' %(max_interval))
    restart_inp_file.write('    md_steps %d\n' %(md_steps))
    if ( os.path.exists(restart_file_name_abs) ):
      restart_inp_file.write('    restart_step %d\n' %(i))
    else:
      restart_inp_file.write('    restart_step %d\n' %(i-1))
    restart_inp_file.write('    cp2k_exe %s\n' %(cp2k_exe))
    restart_inp_file.write('  &end mix_force_eval\n')
    restart_inp_file.write('&end thermo_int\n')

    restart_inp_file.close()

    target_value = 0.0+(i-1)*incre_value
    if ( os.path.exists(restart_file_name_abs) ):
      upper_file_name_abs = file_tools.upper_file(restart_file_name_abs, work_dir)
      line_num = file_tools.grep_line_num("'VALUES '", upper_file_name_abs, work_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no VALUES keyword in cp2k mixing force_eval file')
        exit()
      else:
        value_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'GROUP_PARTITION'", upper_file_name_abs, work_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no GROUP_PARTITION keyword in cp2k mixing force_eval file')
        exit()
      else:
        group_line_num = line_num[0]
        line = linecache.getline(upper_file_name_abs, group_line_num)
        linecache.clearcache()
        line_split = data_op.split_str(line, ' ', '\n')
        if ( len(line_split) == 3 and all(data_op.eval_str(i) == 1 for i in line_split[1:3]) ):
          proc_num = sum([int(x) for x in line_split[1:3]])
        else:
          log_info.log_error('Input error: GROUP_PARTITION should be two integer in cp2k mixing force_eval file')
          exit()
      cmd = "sed -i '%ds/.*/%s %f/' %s" %(value_line_num, '     VALUES ', target_value, restart_file_name)
      call.call_simple_shell(work_dir, cmd)
      cmd = "cp %s %s" %(restart_file_name, inp_trans_file_name)
      call.call_simple_shell(work_dir, cmd)
      cmd = "cp %s %s" %(restart_file_name, 'RESTART_BAK_FILE')
      call.call_simple_shell(work_dir, cmd)
      cmd = "rm %s %s" %(restart_file_name, upper_file_name_abs)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('cp2k running error for step %d' %(i-1))
      exit()
    cmd = "mpirun -np %d %s %s 1> cp2k.out 2> cp2k.err" %(proc_num, cp2k_exe, inp_trans_file_name)
    subprocess.run(cmd, cwd=work_dir, shell=True)
    cmd = "rm %s" %(inp_trans_file_name)
    call.call_simple_shell(work_dir, cmd)

def mix_force_eval_run(mix_param, work_dir):

  '''
  mix_force_eval_run: kernel function to run mix force_eval

  Args:
    mix_param: dictionary
      mix_param contains parameters for mix force_eval
    work_dir: string
      work_dir is the workding directory for CP2K_kit
  Returns:
    none
  '''

  cp2k_exe = call.call_returns_shell(work_dir, 'which cp2k.popt')
  if ( len(cp2k_exe) == 0 or 'no cp2k.popt in' in cp2k_exe[0] ):
    log_info.log_error('Envrionment error: can not find cp2k executable file, please set the environment for cp2k')
    exit()

  mix_param = check_thermo_int.check_mix_force_eval_inp(mix_param)
  mix_inp_file = mix_param['mix_inp_file']
  max_interval = mix_param['max_interval']
  md_steps = mix_param['md_steps']
  restart_step = mix_param['restart_step']

  print ('MIX_FORCE_EVAL'.center(80, '*'), flush=True)
  print ('Run slow-growth mixing force eval claculation', flush=True)
  print ('%d tasks from initial value %f to endding value %f' %(max_interval, 0.0, 1.0), flush=True)
  mix_force_eval(mix_inp_file, max_interval, md_steps, restart_step, cp2k_exe[0], work_dir)
