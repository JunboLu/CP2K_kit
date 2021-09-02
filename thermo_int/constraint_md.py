#! /usr/env/bin python

import os
import linecache
import subprocess
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import call
from CP2K_kit.tools import file_tools
from CP2K_kit.thermo_int import check_thermo_int

def revise_inp_file(cp2k_inp_file, init_value, micro_steps, colvar_id, work_dir):

  '''
  revise_inp_file: revise cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    init_value: float
      init_value is the initial value.
    micro_steps: int
      micro_step is the steps for each macro step.
    colvar_id: int
      colvar_id is the id of colvar.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    proj_name: string
      proj_name is the name of project in cp2k input file.
    target_str: string
      target_str is the target line exclude the target value.
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  #Get project name
  line_num = file_tools.grep_line_num("'PROJECT'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    proj_line_num = line_num[0]
  else:
    line_num = file_tools.grep_line_num("'PROJECT_NAME'", upper_file_name_abs, work_dir)
    if ( line_num != 0 ):
      proj_line_num = line_num[0]
    else:
      log_info.log_error('Input error: no project name in %s file.' %(cp2k_inp_file))
      exit()

  proj_line = linecache.getline(cp2k_inp_file, proj_line_num)
  proj_line_split = data_op.str_split(proj_line, ' ')
  if ( proj_line_split[len(proj_line_split)-1] == '\n' ):
    proj_line_split.remove(proj_line_split[len(proj_line_split)-1])
  else:
    proj_line_split[len(proj_line_split)-1] = proj_line_split[len(proj_line_split)-1].strip('\n')
  proj_name = proj_line_split[len(proj_line_split)-1]

  #Revise inital target value
  line_num = file_tools.grep_line_num("'TARGET'", upper_file_name_abs, work_dir)
  if ( line_num == 0 ):
    log_info.log_error('Input error: no target in %s file' %(cp2k_inp_file))
    exit()
  else:
    target_line_num = line_num[colvar_id-1]
  line = linecache.getline(cp2k_inp_file, target_line_num)
  line_split = data_op.str_split(line, ' ')
  if ( line_split[len(line_split)-1] == '\n' ):
    line_split.remove(line_split[len(line_split)-1])
  target_str = data_op.comb_list_2_str(line_split[0:(len(line_split)-1)], ' ')

  cmd = "sed -ie '%ds/.*/      %s %f/' %s" %(target_line_num, target_str, init_value, cp2k_inp_file)
  call.call_simple_shell(work_dir, cmd)

  #Revise md steps
  line_num = file_tools.grep_line_num("'STEPS'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    step_line_num = line_num[0]
  else:
    log_info.log_error('Input error: no md steps in %s file' %(cp2k_inp_file))

  cmd = "sed -ie '%ds/.*/     STEPS %d/' %s" %(step_line_num, micro_steps, cp2k_inp_file)
  call.call_simple_shell(work_dir, cmd)

  #Revise basis set file name
  line_num = file_tools.grep_line_num("'BASIS_SET_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    basis_line_num = line_num[0]
    basis_line = linecache.getline(cp2k_inp_file, basis_line_num)
    basis_line_split = data_op.str_split(basis_line, ' ')
    if ( basis_line_split[len(basis_line_split)-1] == '\n' ):
      basis_line_split.remove(basis_line_split[len(basis_line_split)-1])
    else:
      basis_line_split[len(basis_line_split)-1] = basis_line_split[len(basis_line_split)-1].strip('\n')
    basis_file_name_abs = os.path.abspath(basis_line_split[len(basis_line_split)-1])
    if ( os.path.exists(basis_file_name_abs) ):
      basis_file_name_abs_split = data_op.str_split(basis_file_name_abs, '/')
      for i in range(len(basis_file_name_abs_split)):
        basis_file_name_abs_split[i] = '\/'+basis_file_name_abs_split[i]
      basis_file_name_abs_trans = data_op.comb_list_2_str(basis_file_name_abs_split, '')
      cmd = "sed -ie '%ds/.*/     BASIS_SET_FILE_NAME %s/' %s" %(basis_line_num, basis_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
  else:
    log_info.log_error('Input error: no basis_set_file_name keyword in %s file' %(cp2k_inp_file))
    exit()

  #Revise potential file name
  line_num = file_tools.grep_line_num("'POTENTIAL_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    pot_line_num = line_num[0]
    pot_line = linecache.getline(cp2k_inp_file, pot_line_num)
    pot_line_split = data_op.str_split(pot_line, ' ')
    if ( pot_line_split[len(pot_line_split)-1] == '\n' ):
      pot_line_split.remove(pot_line_split[len(pot_line_split)-1])
    else:
      pot_line_split[len(pot_line_split)-1] = pot_line_split[len(pot_line_split)-1].strip('\n')
    pot_file_name_abs = os.path.abspath(pot_line_split[len(pot_line_split)-1])
    if ( os.path.exists(pot_file_name_abs) ):
      pot_file_name_abs_split = data_op.str_split(pot_file_name_abs, '/')
      for i in range(len(pot_file_name_abs_split)):
        pot_file_name_abs_split[i] = '\/'+pot_file_name_abs_split[i]
      pot_file_name_abs_trans = data_op.comb_list_2_str(pot_file_name_abs_split, '')
      cmd = "sed -ie '%ds/.*/     BASIS_SET_FILE_NAME %s/' %s" %(pot_line_num, pot_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)

  #Revise coordination file name
  line_num = file_tools.grep_line_num("'COORD_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    coord_line_num = line_num[0]
    coord_line = linecache.getline(cp2k_inp_file, coord_line_num)
    coord_line_split = data_op.str_split(coord_line, ' ')
    if ( coord_line_split[len(coord_line_split)-1] == '\n' ):
      coord_line_split.remove(coord_line_split[len(coord_line_split)-1])
    else:
      coord_line_split[len(coord_line_split)-1] = coord_line_split[len(coord_line_split)-1].strip('\n')
    coord_file_name_abs = os.path.abspath(coord_line_split[len(coord_line_split)-1])
    if ( os.path.exists(coord_file_name_abs) ):
      coord_file_name_abs_split = data_op.str_split(coord_file_name_abs, '/')
      for i in range(len(coord_file_name_abs_split)):
        coord_file_name_abs_split[i] = '\/'+coord_file_name_abs_split[i]
      coord_file_name_abs_trans = data_op.comb_list_2_str(coord_file_name_abs_split, '')
      cmd = "sed -ie '%ds/.*/     COORD_FILE_NAME %s/' %s" %(coord_line_num, coord_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('%s in cp2k input file does not exist' %(coord_line_split[len(coord_line_split)-1]))
      exit()

  #Revise include file name
  line_num = file_tools.grep_line_num("'@INCLUDE'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    for i in range(len(line_num)):
      include_line_num = line_num[i]
      include_line = linecache.getline(cp2k_inp_file, include_line_num)
      include_line_split = data_op.str_split(include_line, ' ')
      if ( include_line_split[len(include_line_split)-1] == '\n' ):
        include_line_split.remove(include_line_split[len(include_line_split)-1])
      else:
        include_line_split[len(include_line_split)-1] = include_line_split[len(include_line_split)-1].strip('\n')
      include_file_name_abs = os.path.abspath(include_line_split[len(include_line_split)-1])
      if ( os.path.exists(include_file_name_abs) ):
        include_file_name_abs_split = data_op.str_split(include_file_name_abs, '/')
        for i in range(len(include_file_name_abs_split)):
          include_file_name_abs_split[i] = '\/'+include_file_name_abs_split[i]
        include_file_name_abs_trans = data_op.comb_list_2_str(include_file_name_abs_split, '')
        cmd = "sed -ie '%ds/.*/     COORD_FILE_NAME %s/' %s" %(include_line_num, include_file_name_abs_trans, cp2k_inp_file)
      else:
        log_info.log_error('%s in cp2k input file does not exist' %(include_line_split[len(include_line_split)-1]))
        exit()

  linecache.clearcache()

  cmd = "rm %s %se*" %(upper_file_name_abs, cp2k_inp_file)
  call.call_simple_shell(work_dir, cmd)

  return proj_name, target_str

def finite_points(cp2k_inp_file, init_value, end_value, macro_steps, micro_steps, \
                  restart_step, colvar_id, cp2k_exe, cp2k_mpi_num, work_dir):

  '''
  finit_points: run finite-points constraind md calculations.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    init_value: float
      init_value is the initial value.
    end_value: float
      end_value is the endding value.
    macro_steps: int
      macro_steps is the totals steps.
    micro_steps: int
      micro_step is the steps for each macro step.
    restart_step: int
      restart_step is the restarting step.
    colvar_id: int
      colvar_id is the id of colvar.
    cp2k_exe: string
      cp2k_exe is the cp2k executable file.
    cp2k_mpi_num: int
      cp2k_mpi_num is the cp2k mpi numbers.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  proj_name, target_str = revise_inp_file(cp2k_inp_file, init_value, micro_steps, colvar_id, work_dir)

  incre_value = float((end_value-init_value)/macro_steps)

  restart_file_name = ''.join((proj_name, '-1.restart'))
  lag_file_name = ''.join((proj_name, '-1.LagrangeMultLog'))

  for i in range(restart_step, macro_steps+2, 1):
    print (''.join(('Task ', str(i))).center(80, '*'), flush=True)
    #Write cp2k_kit restart file
    while True:
      restart_inp_file_name_abs = ''.join((work_dir, '/CP2K_KIT.restart'))
      if ( os.path.exists(restart_inp_file_name_abs) ):
        cmd = 'rm %s' %(restart_inp_file_name_abs)
        call.call_simple_shell(work_dir, cmd)
      else:
        break
    restart_inp_file = open(restart_inp_file_name_abs, 'w')
    restart_inp_file.write('&global\n')
    restart_inp_file.write('  run_type thermo_int\n')
    restart_inp_file.write('  thermo_type constraint_md\n')
    restart_inp_file.write('&end global\n')
    restart_inp_file.write('\n')

    restart_inp_file.write('&thermo_int\n')
    restart_inp_file.write('  &constraint_md\n')
    restart_inp_file.write('    run_type finite_points\n')
    restart_inp_file.write('    cp2k_inp_file %s\n' %(cp2k_inp_file))
    restart_inp_file.write('    init_value %f\n' %(init_value))
    restart_inp_file.write('    end_value %f\n' %(end_value))
    restart_inp_file.write('    macro_steps %d\n' %(macro_steps))
    restart_inp_file.write('    micro_steps %d\n' %(micro_steps))
    restart_inp_file.write('    restart_step %d\n' %(i))
    restart_inp_file.write('    cp2k_exe %s\n' %(cp2k_exe))
    restart_inp_file.write('    cp2k_mpi_num %d\n' %(cp2k_mpi_num))
    restart_inp_file.write('  &end constraint_md\n')
    restart_inp_file.write('&end thermo_int\n')

    restart_inp_file.close()

    target_value = init_value+(i-1)*incre_value
    task_dir = ''.join((work_dir, '/task_', str(i)))
    restart_file_name_abs = ''.join((task_dir, '/', restart_file_name))
    lag_file_name_abs = ''.join((task_dir, '/', lag_file_name))
    if ( os.path.exists(task_dir) and os.path.exists(restart_file_name_abs) ):
      cmd_1 = "grep %s %s" %("'STEP_START_VAL'", restart_file_name_abs)
      line_1 = call.call_returns_shell(task_dir, cmd_1)
      line_1_split = data_op.str_split(line_1[0], ' ')
      if ( line_1_split[len(line_1_split)-1] == '\n' ):
        line_1_split.remove(line_1_split[len(line_1_split)-1])
      else:
        line_1_split[len(line_1_split)-1] = line_1_split[len(line_1_split)-1].strip('\n')
      cmd_2 = "grep %s %s" %("'STEP_START_VAL'", 'input.inp')
      line_2 = call.call_returns_shell(task_dir, cmd_2)
      if ( len(line_2) == 0 ):
        num = int(line_1_split[len(line_1_split)-1])
      else:
        line_2_split = data_op.str_split(line_2[0], ' ')
        if ( line_2_split[len(line_2_split)-2] == '\n' ):
          line_2_split.remove(line_2_split[len(line_2_split)-1])
        else:
          line_2_split[len(line_2_split)-1] = line_2_split[len(line_2_split)-1].strip('\n')
        num = int(line_1_split[len(line_1_split)-1]) - int(line_2_split[len(line_2_split)-1])
      whole_line_num = len(open(lag_file_name_abs, 'r').readlines())
      if ( whole_line_num > 2*num ):
        cmd = "sed -ie '%d,%dd' %s" %(2*num+1, whole_line_num, lag_file_name)
        call.call_simple_shell(task_dir, cmd)
        cmd = "rm %s" %(''.join((lag_file_name, 'e')))
        call.call_simple_shell(task_dir, cmd)
      remain_num = micro_steps-num
      step_line_num = file_tools.grep_line_num("'STEPS'", restart_file_name_abs, work_dir)[0]
      cmd = "sed -ie '%ds/.*/     STEPS %d/' %s " %(step_line_num, remain_num, restart_file_name)
      call.call_simple_shell(task_dir, cmd)
      cmd = "cp %s %s" %(''.join((proj_name, '-RESTART.wfn.bak-1')), ''.join((proj_name, '-RESTART.wfn.bak')))
      call.call_simple_shell(task_dir, cmd)
      print ('Success: generate input files for task %d' %(i), flush=True)
      cmd = "mpirun -np %d %s %s 1>> cp2k.out 2>> cp2k.err" %(cp2k_mpi_num, cp2k_exe, restart_file_name)
      subprocess.run(cmd, cwd=task_dir, shell=True)
      label_line_num = file_tools.grep_line_num("'PROGRAM ENDED AT'", 'cp2k.out', task_dir)[0]
      if ( label_line_num == 0 ):
        log_info.log_error('cp2k running error at %d step' %(i))
        exit()
      else:
        print ('Success: constraint molecular dynamics for task %d' %(i), flush=True)
    else:
      if ( not os.path.exists(task_dir) ):
        cmd = "mkdir %s" %(task_dir)
        call.call_simple_shell(work_dir, cmd)
      else:
        cmd = "rm -rf %s" %(task_dir)
        call.call_simple_shell(work_dir, cmd)
        cmd = "mkdir %s" %(task_dir)
        call.call_simple_shell(work_dir, cmd)
      if ( i == 1 ):
        cmd = "cp %s %s" %(cp2k_inp_file, task_dir)
        call.call_simple_shell(work_dir, cmd)
        cmd = "cp %s %s" %(cp2k_inp_file, 'input.inp')
        call.call_simple_shell(task_dir, cmd)
        print ('Success: generate input files for task %d' %(i), flush=True)
        cmd = "mpirun -np %d %s %s 1> cp2k.out 2> cp2k.err" %(cp2k_mpi_num, cp2k_exe, cp2k_inp_file)
        subprocess.run(cmd, cwd=task_dir, shell=True)
        label_line_num = file_tools.grep_line_num("'PROGRAM ENDED AT'", 'cp2k.out', task_dir)
        if ( label_line_num == 0 ):
          log_info.log_error('cp2k running error at %d step' %(i))
          exit()
        else:
          print ('Success: constraint molecular dynamics for task %d' %(i), flush=True)
      else:
        prev_dir = ''.join((work_dir, '/task_', str(i-1)))
        cmd = "cp %s %s %s" %(restart_file_name, ''.join((proj_name, '-RESTART.wfn')), task_dir)
        call.call_simple_shell(prev_dir, cmd)
        target_line_num = file_tools.grep_line_num("'TARGET'", restart_file_name, task_dir)[colvar_id-1]
        cmd = "sed -ie '%ds/.*/       %s %f/' %s" %(target_line_num, target_str, target_value, restart_file_name)
        call.call_simple_shell(task_dir, cmd)
        step_line_num = file_tools.grep_line_num("'STEPS'", restart_file_name, task_dir)[0]
        cmd = "sed -ie '%ds/.*/     STEPS %d/' %s " %(step_line_num, micro_steps, restart_file_name)
        call.call_simple_shell(task_dir, cmd)
        cmd = "rm %s*" %(''.join((restart_file_name, 'e')))
        call.call_simple_shell(task_dir, cmd)
        cmd = "cp %s %s" %(restart_file_name, 'input.inp')
        call.call_simple_shell(task_dir, cmd)
        print ('Success: generate input files for task %d' %(i), flush=True)
        cmd = "mpirun -np %d %s %s 1>> cp2k.out 2>> cp2k.err" %(cp2k_mpi_num, cp2k_exe, restart_file_name)
        subprocess.run(cmd, cwd=task_dir, shell=True)
        label_line_num = file_tools.grep_line_num("'PROGRAM ENDED AT'", 'cp2k.out', task_dir)
        if ( label_line_num == 0 ):
          log_info.log_error('cp2k running error at %d step' %(i))
          exit()
        else:
          print ('Success: constraint molecular dynamics for task %d' %(i), flush=True)

def slow_growth(cp2k_inp_file, init_value, end_value, macro_steps, micro_steps, \
                restart_step, colvar_id, cp2k_exe, cp2k_mpi_num, work_dir):

  '''
  slow_growth: run slow growth constraind md calculations.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    init_value: float
      init_value is the initial value.
    end_value: float
      end_value is the endding value.
    macro_steps: int
      macro_steps is the totals steps.
    micro_steps: int
      micro_step is the steps for each macro step.
    restart_step: int
      restart_step is the restarting step.
    colvar_id: int
      colvar_id is the id of colvar.
    cp2k_exe: string
      cp2k_exe is the cp2k executable file.
    cp2k_mpi_num: int
      cp2k_mpi_num is the cp2k mpi numbers.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  proj_name, target_str = revise_inp_file(cp2k_inp_file, init_value, micro_steps, colvar_id, work_dir)

  incre_value = float((end_value-init_value)/macro_steps)

  restart_file_name = ''.join((proj_name, '-1.restart'))
  restart_file_name_abs = ''.join((work_dir, '/', restart_file_name))
  lag_file_name = ''.join((proj_name, '-1.LagrangeMultLog'))
  lag_file_name_abs = ''.join((work_dir, '/', lag_file_name))

  if ( restart_step == 1 ):
    cmd = "cp %s %s" %(cp2k_inp_file, restart_file_name)
    call.call_simple_shell(work_dir, cmd)

  m = 0
  while True:
    inp_trans_file_name = ''.join(('input_cmd', str(m), '.inp'))
    inp_trans_file_name_abs = ''.join((work_dir, '/', inp_trans_file_name))
    if ( os.path.exists(inp_trans_file_name_abs) ):
      m = m+1
    else:
      break

  for i in range(restart_step, macro_steps+2, 1):
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

      whole_line_num = len(open(lag_file_name_abs, 'r').readlines())
      if ( whole_line_num > 2*(restart_step-1) ):
        cmd = "sed -ie '%d,%dd' %s" %(2*(restart_step-1)+1, whole_line_num, lag_file_name)
        call.call_simple_shell(work_dir, cmd)
        cmd = "rm %s" %(''.join((lag_file_name, 'e')))
        call.call_simple_shell(work_dir, cmd)

    restart_inp_file = open(restart_inp_file_name_abs, 'w')
    restart_inp_file.write('&global\n')
    restart_inp_file.write('  run_type thermo_int\n')
    restart_inp_file.write('  thermo_type constraint_md\n')
    restart_inp_file.write('&end global\n')
    restart_inp_file.write('\n')

    restart_inp_file.write('&thermo_int\n')
    restart_inp_file.write('  &constraint_md\n')
    restart_inp_file.write('    run_type slow_growth\n')
    restart_inp_file.write('    cp2k_inp_file %s\n' %(cp2k_inp_file))
    restart_inp_file.write('    init_value %f\n' %(init_value))
    restart_inp_file.write('    end_value %f\n' %(end_value))
    restart_inp_file.write('    macro_steps %d\n' %(macro_steps))
    restart_inp_file.write('    micro_steps %d\n' %(micro_steps))
    if ( os.path.exists(restart_file_name_abs) ):
      restart_inp_file.write('    restart_step %d\n' %(i))
    else:
      restart_inp_file.write('    restart_step %d\n' %(i-1))
    restart_inp_file.write('    cp2k_exe %s\n' %(cp2k_exe))
    restart_inp_file.write('    cp2k_mpi_num %d\n' %(cp2k_mpi_num))
    restart_inp_file.write('  &end constraint_md\n')
    restart_inp_file.write('&end thermo_int\n')

    restart_inp_file.close()

    target_value = init_value+(i-1)*incre_value

    if ( os.path.exists(restart_file_name_abs) ):
      upper_file_name_abs = file_tools.upper_file(restart_file_name_abs, work_dir)
      target_line_num = file_tools.grep_line_num("'TARGET'", upper_file_name_abs, work_dir)[colvar_id-1]
      cmd = "sed -ie '%ds/.*/      %s %f/' %s" %(target_line_num, target_str, target_value, restart_file_name_abs)
      call.call_simple_shell(work_dir, cmd)
      cmd = "cp %s %s" %(restart_file_name, inp_trans_file_name)
      call.call_simple_shell(work_dir, cmd)
      cmd = "cp %s %s" %(restart_file_name, 'RESTART_BAK_FILE')
      call.call_simple_shell(work_dir, cmd)
      cmd = "rm %s* %s" %(restart_file_name, upper_file_name_abs)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('cp2k running error in step %d' %(i-1))
      exit()
    cmd = "mpirun -np %d %s %s 1> cp2k.out 2> cp2k.err" %(cp2k_mpi_num, cp2k_exe, inp_trans_file_name)
    subprocess.run(cmd, cwd=work_dir, shell=True)
    cmd = "rm %s" %(inp_trans_file_name)
    call.call_simple_shell(work_dir, cmd)

def constraint_md_run(cmd_param, work_dir):

  '''
  cmd_run: kernel function to run constraint md.

  Args:
    cmd_param: dictionary
      cmd_param contains parameters for constraint md.
    work_dir: string
      work_dir is the workding directory of CP2K_kit.
  Returns:
    none
  '''

  cmd_param = check_thermo_int.check_constraint_md_inp(cmd_param)

  run_type = cmd_param['run_type']
  cp2k_inp_file = cmd_param['cp2k_inp_file']
  init_value = cmd_param['init_value']
  end_value = cmd_param['end_value']
  macro_steps = cmd_param['macro_steps']
  micro_steps = cmd_param['micro_steps']
  restart_step = cmd_param['restart_step']
  cp2k_exe = cmd_param['cp2k_exe']
  cp2k_mpi_num = cmd_param['cp2k_mpi_num']
  colvar_id = cmd_param['colvar_id']

  print ('CONSTRAINT_MD'.center(80, '*'), flush=True)

  if ( run_type == 'slow_growth' ):
    print ('Run slow-growth constraint molecular dynamics', flush=True)
    print ('%d tasks from initial value %f to endding value %f' %(macro_steps+1, init_value, end_value), flush=True)
    slow_growth(cp2k_inp_file, init_value, end_value, macro_steps, micro_steps, \
                restart_step, colvar_id, cp2k_exe, cp2k_mpi_num, work_dir)

  elif ( run_type == 'finite_points' ):
    print ('Run finite-points constraint molecular dynamics', flush=True)
    print ('%d tasks from initial value %f to endding value %f\n' %(macro_steps+1, init_value, end_value), flush=True)
    finite_points(cp2k_inp_file, init_value, end_value, macro_steps, micro_steps, \
                  restart_step, colvar_id, cp2k_exe, cp2k_mpi_num, work_dir)
