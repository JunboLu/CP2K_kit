#! /usr/env/bin python

import os
import copy
import math
import linecache
import numpy as np
from CP2K_kit.tools import *

def gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, sys_id, coord, box, train_stress):

  '''
  gen_cp2kfrc_file: generate cp2k input file for force calculation

  Args:
    cp2k_param: dictionary
      cp2k_param includes parameters of cp2k force calculation.
    work_dir: string
      work_dir is the workding directory of CP2K_kit.
    iter_id: int
      iter_id is current iteration number.
    sys_id: int
      sys_id is id of system.
    coord: 2-d string list, dim = (num of atoms) * 4
      example: [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
                 [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
    box: 3-d string list
      example: [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
                 [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
    train_stress: bool
      train_stress is whether we need to get stress tensor
  Returns:
    none
  '''

  iter_dir = ''.join((work_dir + '/', 'iter_', str(iter_id)))
  cp2k_calc_dir = ''.join((iter_dir, '/03.cp2k_calc'))

  cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_' + str(sys_id)))
  if ( not os.path.exists(cp2k_sys_dir) ):
    cmd = "mkdir %s" % (''.join(('sys_', str(sys_id))))
    call.call_simple_shell(cp2k_calc_dir, cmd)

  atoms = []
  if ( len(coord) != 0 ):
    std_inp_file_name_abs = ''.join((cp2k_sys_dir, '/cp2k_std.inp'))
    std_inp_file = open(std_inp_file_name_abs, 'w')

    std_inp_file.write('&Global\n')
    std_inp_file.write('  PRINT_LEVEL LOW\n')
    std_inp_file.write('  PROJECT_NAME cp2k\n')
    std_inp_file.write('  RUN_TYPE ENERGY_FORCE\n')
    std_inp_file.write('&END Global\n')
    std_inp_file.write('\n')

    std_inp_file.write('&FORCE_EVAL\n')
    std_inp_file.write('  METHOD Quickstep\n')
    std_inp_file.write('  STRESS_TENSOR ANALYTICAL\n')

    cp2k_inp_file = cp2k_param['cp2k_inp_file'][sys_id]
    if ( cp2k_inp_file != 'none' ):
      cp2k_inp_bak_file_name_abs = ''.join((cp2k_sys_dir, '/cp2k_bak.inp'))
      call.call_simple_shell(work_dir, "cp %s %s" %(cp2k_inp_file, cp2k_inp_bak_file_name_abs))
      revise_cp2k_inp.revise_basis_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_dir)
      revise_cp2k_inp.revise_pot_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_dir)
      revise_cp2k_inp.revise_dftd3_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_dir)
      revise_cp2k_inp.revise_rvv10_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_dir)

      cp2k_inp_file_upper = file_tools.upper_file(cp2k_inp_bak_file_name_abs, cp2k_sys_dir)
      cp2k_inp_file_space = file_tools.space_file(cp2k_inp_bak_file_name_abs, ' ', cp2k_sys_dir)

      line_num = file_tools.grep_line_num("'&DFT'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &DFT keywords in %s file' %(cp2k_inp_file))
        exit()
      else:
        dft_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'&END DFT'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &END DFT keywords in %s file' %(cp2k_inp_file))
        exit()
      else:
        end_dft_line_num = line_num[0]
      for i in range(end_dft_line_num-dft_line_num+1):
        line = linecache.getline(cp2k_inp_bak_file_name_abs, dft_line_num+i)
        std_inp_file.write(line)
    else:
      std_inp_file.write('  &DFT\n')

      basis_file_name = cp2k_param['basis_set_file_name']
      psp_file_name = cp2k_param['potential_file_name']
      charge = cp2k_param['charge']
      spin = cp2k_param['multiplicity']
      cutoff = cp2k_param['cutoff']
      std_inp_file.write(''.join(('    BASIS_SET_FILE_NAME ', basis_file_name, '\n')))
      std_inp_file.write(''.join(('    POTENTIAL_FILE_NAME ', psp_file_name, '\n')))
      std_inp_file.write(''.join(('    CHARGE ', charge, '\n')))
      std_inp_file.write(''.join(('    MULTIPLICITY ', spin, '\n')))
      std_inp_file.write('    &MGRID\n')
      std_inp_file.write(''.join(('      CUTOFF ', cutoff, '\n')))
      std_inp_file.write('      NGRIDS 4\n')
      std_inp_file.write('      REL_CUTOFF 50\n')
      std_inp_file.write('    &END MGRID\n')
      std_inp_file.write('    &QS\n')
      std_inp_file.write('      METHOD GPW\n')
      std_inp_file.write('      EXTRAPOLATION ASPC\n')
      std_inp_file.write('    &END QS\n')
      poisson_periodic = cp2k_param['poisson_periodic']
      std_inp_file.write('    &POISSON\n')
      if ( poisson_periodic == 'NONE' ):
        std_inp_file.write('      POISSON_SOLVER WAVELET\n')
      std_inp_file.write('      PERIODIC %s\n' %(poisson_periodic))
      std_inp_file.write('    &END POISSON\n')
      std_inp_file.write('    &SCF\n')
      std_inp_file.write('      MAX_SCF 100\n')
      std_inp_file.write('      SCF_GUESS RESTART\n')
      std_inp_file.write('      EPS_SCF 1.0E-6\n')
      std_inp_file.write('      CHOLESKY INVERSE_DBCSR\n')
      std_inp_file.write('      &OUTER_SCF\n')
      std_inp_file.write('        MAX_SCF 6\n')
      std_inp_file.write('        EPS_SCF 1.0E-6\n')
      std_inp_file.write('      &END OUTER_SCF\n')
      std_inp_file.write('      &OT\n')
      std_inp_file.write('        MINIMIZER CG\n')
      std_inp_file.write('        PRECONDITIONER FULL_ALL\n')
      std_inp_file.write('      &END OT\n')
      std_inp_file.write('    &END SCF\n')

      std_inp_file.write('    &XC\n')
      xc_func = cp2k_param['xc_functional']
      std_inp_file.write(''.join(('      &XC_FUNCTIONAL ', xc_func, '\n')))
      std_inp_file.write('      &END XC_FUNCTIONAL\n')
      if cp2k_param['dftd3'] :
        a_vec = np.array([float(x) for x in box[i][0]])
        b_vec = np.array([float(x) for x in box[i][1]])
        c_vec = np.array([float(x) for x in box[i][2]])
        a, b, c, alpha, beta, gamma = get_cell.get_cell_const(a_vec, b_vec, c_vec)
        l = max(a,b,c)
        d3_file_name = cp2k_param['dftd3_file']
        std_inp_file.write('      &VDW_POTENTIAL\n')
        std_inp_file.write('        POTENTIAL_TYPE PAIR_POTENTIAL\n')
        std_inp_file.write('        &PAIR_POTENTIAL\n')
        std_inp_file.write(''.join(('          R_CUTOFF ', str(l/2.0), '\n')))
        std_inp_file.write('          TYPE DFTD3\n')
        std_inp_file.write(''.join(('          PARAMETER_FILE_NAME ', d3_file_name, '\n')))
        std_inp_file.write(''.join(('          REFERENCE_FUNCTIONAL ', xc_func, '\n')))
        std_inp_file.write('        &END PAIR_POTENTIAL\n')
        std_inp_file.write('      &END VDW_POTENTIAL\n')
      std_inp_file.write('    &END XC\n')
      std_inp_file.write('  &END DFT\n')

    std_inp_file.write('  &PRINT\n')
    std_inp_file.write('    &FORCES\n')
    std_inp_file.write('      FILENAME\n')
    std_inp_file.write('    &END FORCES\n')
    if train_stress:
      std_inp_file.write('    &STRESS_TENSOR\n')
      std_inp_file.write('      FILENAME\n')
      std_inp_file.write('    &END STRESS_TENSOR\n')
    std_inp_file.write('  &END PRINT\n')

    std_inp_file.write('  &SUBSYS\n')
    std_inp_file.write('    &CELL\n')
    std_inp_file.write('      @include box\n')

    if ( cp2k_inp_file != 'none' ):
      line_num = file_tools.grep_line_num("'&CELL'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &CELL keyword in cp2k input file')
        exit()
      else:
        cell_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'&END CELL'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &END CELL keyword in cp2k input file')
        exit()
      else:
        end_cell_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'PERIODIC'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( line_num == 0 ):
        std_inp_file.write('      PERIODIC XYZ\n')
      else:
        for i in line_num:
          if ( i > cell_line_num and i < end_cell_line_num ):
            cell_periodic_line_num = i
        if 'cell_periodic_line_num' in locals():
          line = linecache.getline(cp2k_inp_bak_file_name_abs, cell_periodic_line_num)
          std_inp_file.write(line)
        else:
          std_inp_file.write('      PERIODIC XYZ\n')
    else:
      cell_periodic = cp2k_param['cell_periodic']
      std_inp_file.write('      PERIODIC %s\n' %(cell_periodic))

    std_inp_file.write('    &END CELL\n')
    std_inp_file.write('    &COORD\n')
    std_inp_file.write('      @include coord\n')
    std_inp_file.write('    &END COORD\n')

    if ( cp2k_inp_file != 'none' ):
      for i in range(len(coord[0])):
        atoms.append(coord[0][i][0])
      coord_atoms_type = data_op.list_replicate(atoms)
      kind_atoms_type = []
      kind_line_num = file_tools.grep_line_num("'&KIND'", cp2k_inp_file_space, cp2k_sys_dir)
      end_kind_line_num = file_tools.grep_line_num("'&END KIND'", cp2k_inp_file_space, cp2k_sys_dir)
      if ( len(kind_line_num) != len(end_kind_line_num) ):
        log_info.log_error('Please use &KIND and &END KIND keywords define each atom kind in cp2k input file')
        eixt()
      else:
        for i in range(len(kind_line_num)):
          kind_line = linecache.getline(cp2k_inp_bak_file_name_abs, kind_line_num[i])
          kind_line_split = data_op.split_str(kind_line, ' ', '\n')
          kind_atoms_type.append(kind_line_split[len(kind_line_split)-1])
          for j in range(end_kind_line_num[i]-kind_line_num[i]+1):
            line = linecache.getline(cp2k_inp_bak_file_name_abs, kind_line_num[i]+j)
            std_inp_file.write(line)
      linecache.clearcache()
      for atoms_type in coord_atoms_type:
        if ( atoms_type not in kind_atoms_type ):
          log_info.log_error('Input error: no kind setting for %s atom in cp2k input file, please check cp2k input file' %(atoms_type))
          exit()

      call.call_simple_shell(cp2k_sys_dir, 'rm %s' %(cp2k_inp_bak_file_name_abs))
      call.call_simple_shell(cp2k_sys_dir, 'rm %s' %(cp2k_inp_file_upper))
      call.call_simple_shell(cp2k_sys_dir, 'rm %s' %(cp2k_inp_file_space))

    else:
      for i in range(len(coord[0])):
        atoms.append(coord[0][i][0])
      atoms_type = data_op.list_replicate(atoms)
      for i in atoms_type:
        std_inp_file.write(''.join(('    &KIND ', i, '\n')))
        if cp2k_param['use_sr_basis']:
          basis_name = ''.join((cp2k_param['basis_level'].upper(), '-MOLOPT-SR-GTH-', atom.get_q_info(i)))
        else:
          basis_name = ''.join((cp2k_param['basis_level'].upper(), '-MOLOPT-GTH-', atom.get_q_info(i)))
        std_inp_file.write(''.join(('      BASIS_SET ', basis_name, '\n')))
        psp_name = ''.join(('GTH-', xc_func, '-', atom.get_q_info(i)))
        std_inp_file.write(''.join(('      POTENTIAL ', psp_name, '\n')))
        std_inp_file.write('    &END KIND\n')

    std_inp_file.write('  &END SUBSYS\n')
    std_inp_file.write('&END FORCE_EVAL')

    std_inp_file.close()

    revise_cp2k_inp.delete_line("'WFN_RESTART_FILE_NAME'", std_inp_file_name_abs, cp2k_sys_dir)
    cp2k_inp_file_upper = file_tools.upper_file(std_inp_file_name_abs, cp2k_sys_dir)
    pot_line_num = file_tools.grep_line_num("'POTENTIAL_FILE_NAME'", cp2k_inp_file_upper, cp2k_sys_dir)[0]
    call.call_simple_shell(cp2k_sys_dir, 'rm %s' %(cp2k_inp_file_upper))
    use_prev_wfn = cp2k_param['use_prev_wfn']

    for i in range(len(coord)):
      cp2k_task_dir = ''.join((cp2k_sys_dir, '/task_', str(i)))
      if ( not os.path.exists(cp2k_task_dir) ):
        cmd = "mkdir %s" %(''.join(('task_', str(i))))
        call.call_simple_shell(cp2k_sys_dir, cmd)

      box_file_name_abs = ''.join((cp2k_task_dir, '/box'))
      box_file = open(box_file_name_abs, 'w')
      a_str = '  '.join((box[i][0][0], box[i][0][1], box[i][0][2]))
      b_str = '  '.join((box[i][1][0], box[i][1][1], box[i][1][2]))
      c_str = '  '.join((box[i][2][0], box[i][2][1], box[i][2][2]))
      box_file.write(''.join(('      A ', a_str, '\n')))
      box_file.write(''.join(('      B ', b_str, '\n')))
      box_file.write(''.join(('      C ', c_str, '\n')))
      box_file.close()

      coord_file_name_abs = ''.join((cp2k_task_dir, '/coord'))
      coord_file = open(coord_file_name_abs, 'w')
      for j in range(len(coord[i])):
        coord_str = '  '.join((coord[i][j][0], coord[i][j][1], coord[i][j][2], coord[i][j][3]))
        coord_file.write(''.join(('    ', coord_str, '\n')))
      coord_file.close()

      call.call_simple_shell(cp2k_task_dir, 'cp %s %s' %(std_inp_file_name_abs, 'input.inp'))
      if use_prev_wfn:
        if ( i != 0 ):
          cmd = "sed -i '%d s/^/    WFN_RESTART_FILE_NAME ..\/task_%d\/cp2k-RESTART.wfn\\n/' input.inp" %(pot_line_num+1, i-1)
          call.call_simple_shell(cp2k_task_dir, cmd)
    call.call_simple_shell(cp2k_task_dir, 'rm %s' %(std_inp_file_name_abs))

def gen_cp2k_task(cp2k_dic, work_dir, iter_id, atoms_type_multi_sys, atoms_num_tot, struct_index, \
                  conv_new_data_num, choose_new_data_num_limit, train_stress, success_ratio, \
                  success_ratio_conv, success_devi_ratio):

  '''
  gen_cp2k_task: generate cp2k tasks based on choosed structure index

  Args:
    cp2k_dic: dictionary
      cp2k_dic includes parameters of cp2k force calculation.
    work_dir: string
      work_dir is the workding directory of CP2K_kit.
    iter_id: int
      iter_id is current iteration number.
    atoms_type_multi_sys: 2-d dictionary, dim = (num of lammps systems) * (num of atom types)
      atoms_type_multi_sys is the atoms type for multi-systems.
      example: {0:{'O':1,'H':2,'N':3},1:{'O':1,'S':2,'N':3}}
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example: {1:192,2:90}
    struct_index: dictionary (both are int)
      example: {0:{0:[2,4,6...], 1:[2,3,4...]}, 1:{0:[3,4,6...], 1:[5,6,7...]}}
                 !                !                    !
               sys_id          task_id              traj_id
    conv_new_data_num: int
      conv_new_data_num is the converge of number of new data.
    choose_new_data_num_limit: int
      choose_new_data_num_limit is the maximum number of choosed new structures.
    train_stress: bool
      train_stress is whether we need to get stress tensor
    success_ratio: float
      success_ratio is the successful ratio for whole systems.
    success_ratio_conv: float
      success_ratio_conv is the converge value of successful ratio.
    success_devi_ratio: float
      success_devi_conv is the deviation successful ratio for whole systems.
  Returns:
    none
  '''

  #copy should be done at first, because the following operators will change it!
  cp2k_param = copy.deepcopy(cp2k_dic)

  iter_dir = ''.join((work_dir + '/', 'iter_', str(iter_id)))
  cp2k_calc_dir = ''.join((iter_dir, '/03.cp2k_calc'))
  if ( not os.path.exists(cp2k_calc_dir) ):
    cmd = "mkdir %s" % ('03.cp2k_calc')
    call.call_simple_shell(iter_dir, cmd)

  lammps_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))
  for key in struct_index:
    box = []
    coord = []

    lammps_sys_dir = ''.join((lammps_dir, '/sys_', str(key)))
    task_num = len(struct_index[key])

    choosed_task = []
    choosed_index_num = []
    for i in range(task_num):
      choosed_index = struct_index[key][i]
      if ( len(choosed_index) < conv_new_data_num and success_ratio >= success_ratio_conv and \
           (success_ratio+success_devi_ratio) >= 0.999 ):
        pass
      else:
        choosed_index_num.append(len(choosed_index))
        choosed_task.append(i)

    choosed_index_num_copy = copy.deepcopy(choosed_index_num)
    if ( sum(choosed_index_num)<choose_new_data_num_limit ):
      pass
    else:
      for i in range(len(choosed_index_num)):
        choosed_index_num[i]=int(choosed_index_num_copy[i]/sum(choosed_index_num_copy)*choose_new_data_num_limit)

    for i in range(len(choosed_task)):
      #data file is from model0 lammps md calculation.
      choosed_index = struct_index[key][choosed_task[i]]
      choosed_index_array = np.array(choosed_index)
      np.random.shuffle(choosed_index_array)
      new_choosed_index = list(choosed_index_array[0:choosed_index_num[i]])

      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(choosed_task[i])))
      if ( len(new_choosed_index) != 0 ):
        sys_task_index = data_op.comb_list_2_str(sorted(new_choosed_index), ' ')
        str_print = 'Choosed index for system %d lammps task %d: %s' %(key, choosed_task[i], sys_task_index)
        str_print = data_op.str_wrap(str_print, 80, '  ')
        print (str_print, flush=True)
        for j in sorted(new_choosed_index):
          box_j = []
          coord_j = []
          data_file_name_abs = ''.join((lammps_sys_task_dir, '/data/data_', str(j), '.lmp'))

          line_4 = linecache.getline(data_file_name_abs, 4)
          line_4_split = data_op.split_str(line_4, ' ')
          line_5 = linecache.getline(data_file_name_abs, 5)
          line_5_split = data_op.split_str(line_5, ' ')
          line_6 = linecache.getline(data_file_name_abs, 6)
          line_6_split = data_op.split_str(line_6, ' ')
          line_7 = linecache.getline(data_file_name_abs, 7)
          line_7_split = data_op.split_str(line_7, ' ')
          Lx = float(line_4_split[1])
          Ly = float(line_5_split[1])
          Lz = float(line_6_split[1])
          xy = float(line_7_split[0])
          xz = float(line_7_split[1])
          yz = float(line_7_split[2])
          Lx_str = numeric.get_as_num_string(Lx)
          Ly_str = numeric.get_as_num_string(Ly)
          Lz_str = numeric.get_as_num_string(Lz)
          xy_str = numeric.get_as_num_string(xy)
          xz_str = numeric.get_as_num_string(xz)
          yz_str = numeric.get_as_num_string(yz)
          zero_str = numeric.get_as_num_string(0.0)
          box_j.append([Lx_str,zero_str,zero_str])
          box_j.append([xy_str,Ly_str,zero_str])
          box_j.append([xz_str,yz_str,Lz_str])
          box.append(box_j)

          atoms_num = atoms_num_tot[key]
          atoms_type_dic = atoms_type_multi_sys[key]
          for k in range(atoms_num):
            line_jk = linecache.getline(data_file_name_abs, k+10+1)
            line_jk_split = data_op.split_str(line_jk, ' ', '\n')
            atom_name = data_op.get_dic_keys(atoms_type_dic, int(line_jk_split[1]))[0]
            x_str = line_jk_split[2]
            y_str = line_jk_split[3]
            z_str = line_jk_split[4]
            coord_j.append([atom_name,x_str,y_str,z_str])

          linecache.clearcache()
          coord.append(coord_j)

    gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, key, coord, box, train_stress)

def run_cp2kfrc(work_dir, iter_id, cp2k_exe, parallel_exe, cp2k_env_file, cp2k_job_per_node, proc_num_per_node, host, ssh, atoms_num_tot):

  '''
  run_force: perform cp2k force calculation

  Args:
    work_dir: string
      work_dir is the workding directory of CP2K_kit.
    iter_id: int
      iter_id is current iteration number.
    cp2k_exe: string
      cp2k_exe is the cp2k executable file.
    parallel_exe: string
      parallel_exe is the parallel executable file.
    cp2k_env_file: string
      cp2k_env_file is the environment setting file of cp2k.
    cp2k_job_per_node: int
      cp2k_job_per_node is the job number of cp2k in each node.
    proc_num_per_node: 1-d int list
      proc_num_per_node is the numbers of processor in each node.
    host: 1-d string list
      host is the name of computational nodes.
    ssh: bool
      ssh is whether we need to ssh.
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example: {1:192,2:90}
  Returns:
    none
  '''

  import subprocess

  cp2k_calc_dir = ''.join((work_dir, '/iter_', str(iter_id), '/03.cp2k_calc'))
  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(cp2k_calc_dir, cmd))

  #check generating cp2k tasks
  check_cp2k_gen = []
  for i in range(sys_num):
   cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
   cmd = "ls | grep %s" %('task_')
   task_num = len(call.call_returns_shell(cp2k_sys_dir, cmd))
   for j in range(task_num):
     cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(j)))
     inp_file_name_abs = ''.join((cp2k_sys_task_dir, '/input.inp'))
     if ( os.path.exists(inp_file_name_abs) and os.path.getsize(inp_file_name_abs) != 0 ):
       check_cp2k_gen.append(0)
     else:
       check_cp2k_gen.append(1)
  if ( all(i == 0 for i in check_cp2k_gen) ):
    str_print = 'Success: generate cp2k tasks in %s' %(cp2k_calc_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating cp2k tasks error, please check iteration %d' %(iter_id))
    exit()

  #run cp2k tasks
  for i in range(sys_num):
    cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" %('task_')
    task_num = len(call.call_returns_shell(cp2k_sys_dir, cmd))
    host_name_proc = []
    for l in range(len(host)):
      host_name_proc.append(''.join((str(proc_num_per_node[l]), '/', host[l])))
    host_info = data_op.comb_list_2_str(host_name_proc, ',')

    calculated_id = 0
    for j in range(task_num):
      frc_file=''.join((cp2k_sys_dir, '/task_', str(j), '/cp2k-1_0.xyz'))
      if ( os.path.exists(frc_file) and len(open(frc_file, 'r').readlines()) == atoms_num_tot[i]+5 ):
        calculated_id = j
      else:
        break

    if ( calculated_id != task_num-1 ):
      run_start = calculated_id
      run_end = run_start+cp2k_job_per_node*len(host)-1
      if ( run_end > task_num-1 ):
        run_end=task_num-1
      cycle = math.ceil((task_num-run_start)/(cp2k_job_per_node*len(host)))

      for j in range(cycle):
        tot_mpi_num_list = []
        for proc_num in proc_num_per_node:
          mpi_num_list = data_op.int_split(proc_num, cp2k_job_per_node)
          for k in range(len(mpi_num_list)):
            if ( mpi_num_list[k]%2 != 0 and mpi_num_list[k]>1 ):
              mpi_num_list[k] = mpi_num_list[k]-1
          tot_mpi_num_list.append(mpi_num_list)
        tot_mpi_num_list = data_op.list_reshape(tot_mpi_num_list)[0:(run_end-run_start+1)]
        mpi_num_str = data_op.comb_list_2_str(tot_mpi_num_list, ' ')
        task_job_list = data_op.gen_list(run_start, run_end, 1)
        task_job_str = data_op.comb_list_2_str(task_job_list, ' ')

        run_1 = '''
#! /bin/bash

task_job="%s"
mpi_num="%s"
direc=%s
parallel_exe=%s

task_job_arr=(${task_job///})
mpi_num_arr=(${mpi_num///})

num=${#task_job_arr[*]}

for ((i=0;i<=num-1;i++));
do
task_job_mpi_num_arr[i]="${task_job_arr[i]} ${mpi_num_arr[i]}"
done
''' %(task_job_str, mpi_num_str, cp2k_sys_dir, parallel_exe)
        if ssh:
          run_2 = '''
for i in "${task_job_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d -S %s --sshdelay 0.1 $direc/produce.sh {} $direc
''' %( cp2k_job_per_node, host_info)
        else:
          run_2 = '''
for i in "${task_job_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d $direc/produce.sh {} $direc
''' %( cp2k_job_per_node)

        produce = '''
#! /bin/bash

source %s

x=$1
direc=$2

x_arr=(${x///})

new_direc=$direc/task_${x_arr[0]}

cd $new_direc
if [ -f "cp2k-1_0.xyz" ]; then
rm cp2k-1_0.xyz
fi
mpirun -np ${x_arr[1]} %s $new_direc/input.inp 1> $new_direc/cp2k.out 2> $new_direc/cp2k.err

converge_info=`grep "SCF run NOT converged" cp2k.out`
if [ $? -eq 0 ]; then
wfn_line=`grep -n "WFN_RESTART_FILE_NAME" input.inp`
if [ $? -eq 0 ]; then
line=`grep -n "WFN_RESTART_FILE_NAME" input.inp | awk -F ":" '{print $1}'`
sed -i ''$line's/.*/    WFN_RESTART_FILE_NAME .\/cp2k-RESTART.wfn/' input.inp
else
line=`grep -n "POTENTIAL_FILE_NAME" input.inp | awk -F ":" '{print $1}'`
sed -i ''$line' s/^/    WFN_RESTART_FILE_NAME .\/cp2k-RESTART.wfn\\n/' input.inp
fi
if [ -f "cp2k-1_0.xyz" ]; then
rm cp2k-1_0.xyz
fi
mpirun -np ${x_arr[1]} %s $new_direc/input.inp 1> $new_direc/cp2k.out 2> $new_direc/cp2k.err
fi
cd %s
''' %(cp2k_env_file, cp2k_exe, cp2k_exe, work_dir)

        run_file_name_abs = ''.join((cp2k_sys_dir, '/run.sh'))
        with open(run_file_name_abs, 'w') as f:
          f.write(run_1+run_2)

        produce_file_name_abs = ''.join((cp2k_sys_dir, '/produce.sh'))
        with open(produce_file_name_abs, 'w') as f:
          f.write(produce)

        subprocess.run('chmod +x run.sh', cwd=cp2k_sys_dir, shell=True)
        subprocess.run('chmod +x produce.sh', cwd=cp2k_sys_dir, shell=True)
        try:
          subprocess.run("bash -c './run.sh'", cwd=cp2k_sys_dir, shell=True)
        except subprocess.CalledProcessError as err:
          log_info.log_error('Running error: %s command running error in %s' %(err.cmd, cp2k_sys_dir))

        run_start = run_start + cp2k_job_per_node*len(host)
        run_end = run_end + cp2k_job_per_node*len(host)
        if ( run_end > task_num-1):
          run_end = task_num-1
    else:
      cmd = "mpirun -np %d %s input.inp 1> cp2k.out 2> cp2k.err" %(sum(proc_num_per_node), cp2k_exe)
      task_dir = ''.join((cp2k_sys_dir, '/task_', str(task_num-1)))
      frc_file_name_abs = ''.join((task_dir, '/cp2k-1_0.xyz'))
      if ( os.path.exists(frc_file_name_abs) ):
        cmd = "rm %s" %(frc_file_name_abs)
        call.call_simple_shell(task_dir, cmd)
      try:
        subprocess.run(cmd, cwd=task_dir, shell=True)
      except subprocess.CalledProcessError as err:
        log_info.log_error('Running error: %s command running error in %s' %(err.cmd, task_dir))

  #check running cp2k tasks
  check_cp2k_run = []
  for i in range(sys_num):
    cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" %('task_')
    task_num = len(call.call_returns_shell(cp2k_sys_dir, cmd))
    check_cp2k_run_i = []
    for j in range(task_num):
      cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(j)))
      frc_file_name_abs = ''.join((cp2k_sys_task_dir, '/cp2k-1_0.xyz'))
      if ( os.path.exists(frc_file_name_abs) and len(open(frc_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+5 ):
        check_cp2k_run_i.append(0)
      else:
        check_cp2k_run_i.append(1)
    check_cp2k_run.append(check_cp2k_run_i)
  if ( all(i == 0 for i in data_op.list_reshape(check_cp2k_run)) ):
    print ('  Success: ab initio force calculations for %d systems by cp2k' %(sys_num), flush=True)
  else:
    for j in range(sys_num):
      failure_task_id = [index for (index,value) in enumerate(check_cp2k_run[j]) if value==1]
      if ( len(failure_task_id) != task_num ):
        failure_task_id_str = data_op.comb_list_2_str(failure_task_id, ' ')
        str_print = '  Warning: ab initio force calculations fail for tasks %s in system %d by cp2k' %(failure_task_id_str, i)
        str_print = data_op.str_wrap(str_print, 80, '  ')
        print (str_print, flush=True)
      else:
        log_info.log_error('Running error: ab initio force calculations running error, please check iteration %d' %(iter_id))
        exit()

if __name__ == '__main__':
  from collections import OrderedDict
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import cp2k_run
  from CP2K_kit.deepff import check_deepff

  work_dir = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md_mtd'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  proc_num = 4
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  #Test gen_cp2kfrc_file function
  coord = [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
            [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
  box = [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
         [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
#  cp2k_run.gen_cp2kfrc_file(cp2k_dic, work_dir, 1, 0, coord, box)

  #Test gen_cp2k_task function
  atoms_type_multi_sys = {0: {'C': 1, 'O': 2}, 1: {'C': 1, 'O': 2}}
  atoms_num_tot = {0:3, 1:3}
  struct_index = OrderedDict([(0, OrderedDict([(0, [237, 264, 275, 291, 331, 367, 422])])), (1, OrderedDict([(0, [])]))])
  conv_new_data_num = 5
  choose_new_data_num_limit = 100
  cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, 17, atoms_type_multi_sys, atoms_num_tot, \
                         struct_index, conv_new_data_num, choose_new_data_num_limit, False)

  #Test run_cp2kfrc function
#  cp2k_run.run_cp2kfrc(work_dir, 0, environ_dic)
