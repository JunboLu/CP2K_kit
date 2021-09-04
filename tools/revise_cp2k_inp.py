#! /usr/env/bin python

import os
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import file_tools
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info

def get_proj_name(cp2k_inp_file, work_dir):

  '''
  get_proj_name: get the project name of cp2k.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
  Returns:
    proj_name: string
      proj_name is the name of project in cp2k input file.
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

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
  linecache.clearcache()
  proj_line_split = data_op.str_split(proj_line, ' ')
  if ( proj_line_split[len(proj_line_split)-1] == '\n' ):
    proj_line_split.remove(proj_line_split[len(proj_line_split)-1])
  else:
    proj_line_split[len(proj_line_split)-1] = proj_line_split[len(proj_line_split)-1].strip('\n')
  proj_name = proj_line_split[len(proj_line_split)-1]

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

  return proj_name

def revise_target_value(cp2k_inp_file, target_value, colvar_id, work_dir):

  '''
  revise_target_value: revise target value in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    target_value: float
      target_value is the target value.
    colvar_id: int
      colvar_id is the id of colvar.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    target_str: string
      target_str is the target line exclude the target value.
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'TARGET'", upper_file_name_abs, work_dir)
  if ( line_num == 0 ):
    log_info.log_error('Input error: no target in %s file' %(cp2k_inp_file))
    exit()
  else:
    target_line_num = line_num[colvar_id-1]
  line = linecache.getline(cp2k_inp_file, target_line_num)
  linecache.clearcache()
  line_split = data_op.str_split(line, ' ')
  if ( line_split[len(line_split)-1] == '\n' ):
    line_split.remove(line_split[len(line_split)-1])
  target_str = data_op.comb_list_2_str(line_split[0:(len(line_split)-1)], ' ')

  cmd = "sed -i '%ds/.*/      %s %f/' %s" %(target_line_num, target_str, target_value, cp2k_inp_file)
  call.call_simple_shell(work_dir, cmd)

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

  return target_str

def revise_md_steps(cp2k_inp_file, md_steps, work_dir):

  '''
  revise_md_steps: revise md steps in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    md_steps: int
      micro_step is the steps for each macro step.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'STEPS'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    step_line_num = line_num[0]
  else:
    log_info.log_error('Input error: no md steps in %s file' %(cp2k_inp_file))

  cmd = "sed -i '%ds/.*/     STEPS %d/' %s" %(step_line_num, md_steps, cp2k_inp_file)
  call.call_simple_shell(work_dir, cmd)

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_basis_file_name(cp2k_inp_file, work_dir):

  '''
  revise_basis_file_name: revise basis file name in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'BASIS_SET_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    basis_line_num = line_num[0]
    basis_line = linecache.getline(cp2k_inp_file, basis_line_num)
    linecache.clearcache()
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
      cmd = "sed -i '%ds/.*/     BASIS_SET_FILE_NAME %s/' %s" %(basis_line_num, basis_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
  else:
    log_info.log_error('Input error: no basis_set_file_name keyword in %s file' %(cp2k_inp_file))
    exit()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_pot_file_name(cp2k_inp_file, work_dir):

  '''
  revise_pot_file_name: revise potential file name in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'POTENTIAL_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    pot_line_num = line_num[0]
    pot_line = linecache.getline(cp2k_inp_file, pot_line_num)
    linecache.clearcache()
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
      cmd = "sed -i '%ds/.*/     POTENTIAL_FILE_NAME %s/' %s" %(pot_line_num, pot_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
  else:
    log_info.log_error('Input error: no potential_file_name keyword in %s file' %(cp2k_inp_file))
    exit()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_coord_file_name(cp2k_inp_file, work_dir):

  '''
  revise_coord_file_name: revise the coordination file in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'COORD_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    coord_line_num = line_num[0]
    coord_line = linecache.getline(cp2k_inp_file, coord_line_num)
    linecache.clearcache()
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
      cmd = "sed -i '%ds/.*/     COORD_FILE_NAME %s/' %s" %(coord_line_num, coord_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('%s in cp2k input file does not exist' %(coord_line_split[len(coord_line_split)-1]))
      exit()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_dftd3_file_name(cp2k_inp_file, work_dir):

  '''
  revise_dftd3_file_name: revise dftd3 file name in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'PARAMETER_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    dftd3_line_num = line_num[0]
    dftd3_line = linecache.getline(cp2k_inp_file, dftd3_line_num)
    linecache.clearcache()
    dftd3_line_split = data_op.str_split(dftd3_line, ' ')
    if ( dftd3_line_split[len(dftd3_line_split)-1] == '\n' ):
      dftd3_line_split.remove(dftd3_line_split[len(dftd3_line_split)-1])
    else:
      dftd3_line_split[len(dftd3_line_split)-1] = dftd3_line_split[len(dftd3_line_split)-1].strip('\n')
    dftd3_file_name_abs = os.path.abspath(dftd3_line_split[len(dftd3_line_split)-1])
    if ( os.path.exists(dftd3_file_name_abs) ):
      dftd3_file_name_abs_split = data_op.str_split(dftd3_file_name_abs, '/')
      for i in range(len(dftd3_file_name_abs_split)):
        dftd3_file_name_abs_split[i] = '\/'+dftd3_file_name_abs_split[i]
      dftd3_file_name_abs_trans = data_op.comb_list_2_str(dftd3_file_name_abs_split, '')
      cmd = "sed -i '%ds/.*/          PARAMETER_FILE_NAME %s/' %s" %(dftd3_line_num, dftd3_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('%s in cp2k input file does not exist' %(dftd3_line_split[len(dftd3_line_split)-1]))
      exit()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_rvv10_file_name(cp2k_inp_file, work_dir):

  '''
  revise_rvv10_file_name: revise rvv10 file name in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num("'KERNEL_FILE_NAME'", upper_file_name_abs, work_dir)
  if ( line_num != 0 ):
    rvv10_line_num = line_num[0]
    rvv10_line = linecache.getline(cp2k_inp_file, rvv10_line_num)
    linecache.clearcache()
    rvv10_line_split = data_op.str_split(rvv10_line, ' ')
    if ( rvv10_line_split[len(rvv10_line_split)-1] == '\n' ):
      rvv10_line_split.remove(rvv10_line_split[len(rvv10_line_split)-1])
    else:
      rvv10_line_split[len(rvv10_line_split)-1] = rvv10_line_split[len(rvv10_line_split)-1].strip('\n')
    rvv10_file_name_abs = os.path.abspath(rvv10_line_split[len(rvv10_line_split)-1])
    if ( os.path.exists(rvv10_file_name_abs) ):
      rvv10_file_name_abs_split = data_op.str_split(rvv10_file_name_abs, '/')
      for i in range(len(rvv10_file_name_abs_split)):
        rvv10_file_name_abs_split[i] = '\/'+rvv10_file_name_abs_split[i]
      rvv10_file_name_abs_trans = data_op.comb_list_2_str(rvv10_file_name_abs_split, '')
      cmd = "sed -i '%ds/.*/          KERNEL_FILE_NAME %s/' %s" %(rvv10_line_num, rvv10_file_name_abs_trans, cp2k_inp_file)
      call.call_simple_shell(work_dir, cmd)
    else:
      log_info.log_error('%s in cp2k input file does not exist' %(rvv10_line_split[len(rvv10_line_split)-1]))
      exit()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def revise_include_file_name(cp2k_inp_file, work_dir):

  '''
  revise_include_file_name: revise including file name in cp2k input file.

  Args:
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

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
        cmd = "sed -i '%ds/.*/     @include %s/' %s" %(include_line_num, include_file_name_abs_trans, cp2k_inp_file)
        call.call_simple_shell(work_dir, cmd)
      else:
        log_info.log_error('%s in cp2k input file does not exist' %(include_line_split[len(include_line_split)-1]))
        exit()

    linecache.clearcache()

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

def delete_line(keyword, cp2k_inp_file, work_dir):

  '''
  revise_include_file_name: revise including file name in cp2k input file.

  Args:
    keyword: string
      keyword is the keyword in deleted line.
    cp2k_inp_file: string
      cp2k_inp_file is the name of cp2k input file.
    work_dir: string
      work_dir is the working directory of CP2K_kit
  Returns:
    none
  '''

  upper_file_name_abs = file_tools.upper_file(cp2k_inp_file, work_dir)

  line_num = file_tools.grep_line_num(keyword, upper_file_name_abs, work_dir)


  if ( line_num != 0 ):
    keyword_line_num = line_num[0]
    cmd = "sed -i '%dd' %s" %(keyword_line_num, cp2k_inp_file)
    call.call_simple_shell(work_dir, cmd)

  cmd = 'rm %s' %(upper_file_name_abs)
  call.call_simple_shell(work_dir, cmd)

