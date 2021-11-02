#! /usr/env/bin python

import os
import copy
import math
import linecache
import numpy as np
from CP2K_kit.tools import *

def gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, sys_id, task_id, coord, box, train_stress):

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
      sys_id is the id of system.
    task_id: int
      task_id is the id of task.
    coord: 3-d string list, dim = (num of frames)*(num of atoms)*4
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

  cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_' + str(task_id)))
  if ( not os.path.exists(cp2k_sys_task_dir) ):
    cmd = "mkdir %s" % (''.join(('task_', str(task_id))))
    call.call_simple_shell(cp2k_sys_dir, cmd)

  atoms = []
  if ( len(coord) != 0 ):
    std_inp_file_name_abs = ''.join((cp2k_sys_task_dir, '/cp2k_std.inp'))
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
      cp2k_inp_bak_file_name_abs = ''.join((cp2k_sys_task_dir, '/cp2k_bak.inp'))
      call.call_simple_shell(work_dir, "cp %s %s" %(cp2k_inp_file, cp2k_inp_bak_file_name_abs))
      revise_cp2k_inp.revise_basis_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_task_dir)
      revise_cp2k_inp.revise_pot_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_task_dir)
      revise_cp2k_inp.revise_dftd3_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_task_dir)
      revise_cp2k_inp.revise_rvv10_file_name(cp2k_inp_bak_file_name_abs, cp2k_sys_task_dir)

      cp2k_inp_file_upper = file_tools.upper_file(cp2k_inp_bak_file_name_abs, cp2k_sys_task_dir)
      cp2k_inp_file_space = file_tools.space_file(cp2k_inp_bak_file_name_abs, ' ', cp2k_sys_task_dir)

      line_num = file_tools.grep_line_num("'&DFT'", cp2k_inp_file_space, cp2k_sys_task_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &DFT keywords in %s file' %(cp2k_inp_file))
        exit()
      else:
        dft_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'&END DFT'", cp2k_inp_file_space, cp2k_sys_task_dir)
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
      line_num = file_tools.grep_line_num("'&CELL'", cp2k_inp_file_space, cp2k_sys_task_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &CELL keyword in cp2k input file')
        exit()
      else:
        cell_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'&END CELL'", cp2k_inp_file_space, cp2k_sys_task_dir)
      if ( line_num == 0 ):
        log_info.log_error('Input error: no &END CELL keyword in cp2k input file')
        exit()
      else:
        end_cell_line_num = line_num[0]
      line_num = file_tools.grep_line_num("'PERIODIC'", cp2k_inp_file_space, cp2k_sys_task_dir)
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
      kind_line_num = file_tools.grep_line_num("'&KIND'", cp2k_inp_file_space, cp2k_sys_task_dir)
      end_kind_line_num = file_tools.grep_line_num("'&END KIND'", cp2k_inp_file_space, cp2k_sys_task_dir)
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

      call.call_simple_shell(cp2k_sys_task_dir, 'rm %s' %(cp2k_inp_bak_file_name_abs))
      call.call_simple_shell(cp2k_sys_task_dir, 'rm %s' %(cp2k_inp_file_upper))
      call.call_simple_shell(cp2k_sys_task_dir, 'rm %s' %(cp2k_inp_file_space))

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

    revise_cp2k_inp.delete_line("'WFN_RESTART_FILE_NAME'", std_inp_file_name_abs, cp2k_sys_task_dir)
    cp2k_inp_file_upper = file_tools.upper_file(std_inp_file_name_abs, cp2k_sys_task_dir)
    pot_line_num = file_tools.grep_line_num("'POTENTIAL_FILE_NAME'", cp2k_inp_file_upper, cp2k_sys_task_dir)[0]
    call.call_simple_shell(cp2k_sys_task_dir, 'rm %s' %(cp2k_inp_file_upper))
    use_prev_wfn = cp2k_param['use_prev_wfn']

    for i in range(len(coord)):
      cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(i)))
      if ( not os.path.exists(cp2k_sys_task_traj_dir) ):
        cmd = "mkdir %s" %(''.join(('traj_', str(i))))
        call.call_simple_shell(cp2k_sys_task_dir, cmd)

      box_file_name_abs = ''.join((cp2k_sys_task_traj_dir, '/box'))
      box_file = open(box_file_name_abs, 'w')
      a_str = '  '.join((box[i][0][0], box[i][0][1], box[i][0][2]))
      b_str = '  '.join((box[i][1][0], box[i][1][1], box[i][1][2]))
      c_str = '  '.join((box[i][2][0], box[i][2][1], box[i][2][2]))
      box_file.write(''.join(('      A ', a_str, '\n')))
      box_file.write(''.join(('      B ', b_str, '\n')))
      box_file.write(''.join(('      C ', c_str, '\n')))
      box_file.close()

      coord_file_name_abs = ''.join((cp2k_sys_task_traj_dir, '/coord'))
      coord_file = open(coord_file_name_abs, 'w')
      for j in range(len(coord[i])):
        coord_str = '  '.join((coord[i][j][0], coord[i][j][1], coord[i][j][2], coord[i][j][3]))
        coord_file.write(''.join(('    ', coord_str, '\n')))
      coord_file.close()

      call.call_simple_shell(cp2k_sys_task_traj_dir, 'cp %s %s' %(std_inp_file_name_abs, 'input.inp'))
      if use_prev_wfn:
        if ( i != 0 ):
          cmd = "sed -i '%d s/^/    WFN_RESTART_FILE_NAME ..\/traj_%d\/cp2k-RESTART.wfn\\n/' input.inp" %(pot_line_num+1, i-1)
          call.call_simple_shell(cp2k_sys_task_traj_dir, cmd)
    call.call_simple_shell(cp2k_sys_task_dir, 'rm %s' %(std_inp_file_name_abs))

def gen_cp2k_task(cp2k_dic, work_dir, iter_id, atoms_type_multi_sys, atoms_num_tot, struct_index, conv_new_data_num, \
                  choose_new_data_num_limit, train_stress, active_type, success_ratio=0.0, success_devi_ratio=0.0):

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
      train_stress is whether we need to get stress tensor.
    active_type: string
      active_type is the type of active learning.
    success_ratio: float
      success_ratio is the successful ratio for whole systems.
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
    lammps_sys_dir = ''.join((lammps_dir, '/sys_', str(key)))
    task_num = len(struct_index[key])

    choosed_task = []
    choosed_index_num = []
    for i in range(task_num):
      if ( active_type == 'dp_test' ):
        choosed_task.append(i)
      elif ( active_type == 'model_devi' ):
        choosed_index = struct_index[key][i]
        if ( len(choosed_index) < conv_new_data_num and success_ratio >= 0.96 and \
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
      box = []
      coord = []
      #data file is from model0 lammps md calculation.
      if ( active_type == 'dp_test' ):
        choosed_index = struct_index[key][choosed_task[i]]
      elif ( active_type == 'model_devi' ):
        choosed_index = struct_index[key][choosed_task[i]]
        choosed_index_array = np.array(choosed_index)
        np.random.shuffle(choosed_index_array)
        choosed_index = list(choosed_index_array[0:choosed_index_num[i]])

      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(choosed_task[i])))
      if ( len(choosed_index) != 0 ):
        if ( active_type == 'model_devi' ):
          sys_task_index = data_op.comb_list_2_str(sorted(choosed_index), ' ')
          str_print = 'Choosed index for system %d lammps task %d: %s' %(key, choosed_task[i], sys_task_index)
          str_print = data_op.str_wrap(str_print, 80, '  ')
          print (str_print, flush=True)
        for j in range(len(choosed_index)):
          index_j = sorted(choosed_index)[j]
          box_j = []
          coord_j = []
          data_file_name_abs = ''.join((lammps_sys_task_dir, '/data/data_', str(index_j), '.lmp'))

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

        gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, key, choosed_task[i], coord, box, train_stress)

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
