#! /usr/env/bin python

import copy
import linecache
import numpy as np
from CP2K_kit.tools import *

def gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, sys_id, coord, box, get_stress=False):

  '''
  gen_cp2kfrc_file : generate cp2k input file for force calculation

  Args :
    cp2k_param : dictionary
      cp2k_param includes parameters of cp2k force calculation.
    work_dir : string
      work_dir is workding directory.
    iter_id : int
      iter_id is current iteration number.
    sys_id : int
      sys_id is id of system.
    box : 3-d string list
      example : [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
                 [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
    coord : 3-d string list
      example : [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
                 [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
    get_stress : bool
      get_stress is whether we need to get stress tensor
  Returns :
    none
  '''

  iter_dir = ''.join((work_dir + '/', 'iter_', str(iter_id)))
  cp2k_calc_dir = ''.join((iter_dir, '/03.cp2k_calc'))

  cmd = "mkdir %s" % (''.join(('sys_', str(sys_id))))
  call.call_simple_shell(cp2k_calc_dir, cmd)
  cp2k_task_dir = ''.join((cp2k_calc_dir, '/sys_' + str(sys_id)))

  atoms = []
  for i in range(len(coord[0])):
    atoms.append(coord[0][i][0])
  atoms_type = list_dic_op.list_replicate(atoms)

  for i in range(len(coord)):
    cmd = "mkdir %s" %(''.join(('task_', str(i))))
    call.call_simple_shell(cp2k_task_dir, cmd)

    task_dir = ''.join((cp2k_task_dir, '/task_', str(i)))

    task_job_inp = ''.join((task_dir, '/input.inp'))
    inp_file = open(task_job_inp,'w')

    inp_file.write('&Global\n')
    inp_file.write('  PRINT_LEVEL LOW\n')
    inp_file.write('  PROJECT_NAME cp2k\n')
    inp_file.write('  RUN_TYPE ENERGY_FORCE\n')
    inp_file.write('&END Global\n')
    inp_file.write('\n')

    inp_file.write('&FORCE_EVAL\n')
    inp_file.write('  METHOD Quickstep\n')
    inp_file.write('  STRESS_TENSOR ANALYTICAL\n')
    inp_file.write('  &DFT\n')

    basis_file = cp2k_param['basis_set_file_name']
    psp_file = cp2k_param['potential_file_name']
    if ( i ==0 ):
      wfn_file = './cp2k-RESTART.wfn'
    else:
      wfn_file = ''.join(('../task_', str(i-1) , '/cp2k-RESTART.wfn'))
    charge = cp2k_param['charge']
    spin = cp2k_param['multiplicity']
    cutoff = cp2k_param['cutoff']
    inp_file.write(''.join(('    BASIS_SET_FILE_NAME ', basis_file, '\n')))
    inp_file.write(''.join(('    POTENTIAL_FILE_NAME ', psp_file, '\n')))
    inp_file.write(''.join(('    WFN_RESTART_FILE_NAME ', wfn_file, '\n')))
    inp_file.write(''.join(('    CHARGE ', charge, '\n')))
    inp_file.write(''.join(('    MULTIPLICITY ', spin, '\n')))
    inp_file.write('    &MGRID\n')
    inp_file.write(''.join(('      CUTOFF ', cutoff, '\n')))
    inp_file.write('      NGRIDS 4\n')
    inp_file.write('      REL_CUTOFF 50\n')
    inp_file.write('    &END MGRID\n')
    inp_file.write('    &QS\n')
    inp_file.write('      METHOD GPW\n')
    inp_file.write('      EXTRAPOLATION ASPC\n')
    inp_file.write('    &END QS\n')
    inp_file.write('    &POISSON\n')
    inp_file.write('      PERIODIC XYZ\n')
    inp_file.write('    &END POISSON\n')
    inp_file.write('    &SCF\n')
    inp_file.write('      MAX_SCF 100\n')
    inp_file.write('      SCF_GUESS RESTART\n')
    inp_file.write('      EPS_SCF 1.0E-6\n')
    inp_file.write('      CHOLESKY INVERSE_DBCSR\n')
    inp_file.write('      &OUTER_SCF\n')
    inp_file.write('        MAX_SCF 6\n')
    inp_file.write('        EPS_SCF 1.0E-6\n')
    inp_file.write('      &END OUTER_SCF\n')
    inp_file.write('      &OT\n')
    inp_file.write('        MINIMIZER CG\n')
    inp_file.write('        PRECONDITIONER FULL_ALL\n')
    inp_file.write('      &END OT\n')
    inp_file.write('    &END SCF\n')

    inp_file.write('    &XC\n')
    xc_func = cp2k_param['xc_functional']
    inp_file.write(''.join(('      &XC_FUNCTIONAL ', xc_func, '\n')))
    inp_file.write('      &END XC_FUNCTIONAL\n')
    if cp2k_param['dftd3'] :
      a_vec = np.array([float(x) for x in box[i][0]])
      b_vec = np.array([float(x) for x in box[i][1]])
      c_vec = np.array([float(x) for x in box[i][2]])
      a, b, c, alpha, beta, gamma = get_cell.get_cell_const(a_vec, b_vec, c_vec)
      l = max(a,b,c)
      d3_file = cp2k_param['dftd3_file']
      inp_file.write('      &VDW_POTENTIAL\n')
      inp_file.write('        POTENTIAL_TYPE PAIR_POTENTIAL\n')
      inp_file.write('        &PAIR_POTENTIAL\n')
      inp_file.write(''.join(('          R_CUTOFF ', str(l/2.0), '\n')))
      inp_file.write('          TYPE DFTD3\n')
      inp_file.write(''.join(('          PARAMETER_FILE_NAME ', d3_file, '\n')))
      inp_file.write(''.join(('          REFERENCE_FUNCTIONAL ', xc_func, '\n')))
      inp_file.write('        &END PAIR_POTENTIAL\n')
      inp_file.write('      &END VDW_POTENTIAL\n')
    inp_file.write('    &END XC\n')
    inp_file.write('  &END DFT\n')

    inp_file.write('  &PRINT\n')
    inp_file.write('    &FORCES\n')
    inp_file.write('      FILENAME\n')
    inp_file.write('    &END FORCES\n')
    if get_stress:
      inp_file.write('    &STRESS_TENSOR\n')
      inp_file.write('      FILENAME\n')
      inp_file.write('    &END STRESS_TENSOR\n')
    inp_file.write('  &END PRINT\n')

    inp_file.write('  &SUBSYS\n')
    inp_file.write('    &CELL\n')
    a_str = '  '.join((box[i][0][0], box[i][0][1], box[i][0][2]))
    b_str = '  '.join((box[i][1][0], box[i][1][1], box[i][1][2]))
    c_str = '  '.join((box[i][2][0], box[i][2][1], box[i][2][2]))
    inp_file.write(''.join(('      A ', a_str, '\n')))
    inp_file.write(''.join(('      B ', b_str, '\n')))
    inp_file.write(''.join(('      C ', c_str, '\n')))
    inp_file.write('      PERIODIC XYZ\n')
    inp_file.write('    &END CELL\n')
    inp_file.write('    &COORD\n')

    #coord include atoms.
    for j in range(len(coord[i])):
      coord_str = '  '.join((coord[i][j][0], coord[i][j][1], coord[i][j][2], coord[i][j][3]))
      inp_file.write(''.join(('    ', coord_str, '\n')))
    inp_file.write('    &END COORD\n')

    for j in atoms_type:
      inp_file.write(''.join(('    &KIND ', j, '\n')))
      basis_name = ''.join((cp2k_param['basis_level'].upper(), '-MOLOPT-GTH-', atom.get_q_info(j)))
      inp_file.write(''.join(('      BASIS_SET ', basis_name, '\n')))
      psp_name = ''.join(('GTH-', xc_func, '-', atom.get_q_info(j)))
      inp_file.write(''.join(('      POTENTIAL ', psp_name, '\n')))
      inp_file.write('    &END KIND\n')

    inp_file.write('  &END SUBSYS\n')
    inp_file.write('&END FORCE_EVAL')

    inp_file.close()

def gen_cp2k_task(cp2k_dic, work_dir, iter_id, atoms_type_dic_tot, atoms_num_tot, \
                  struct_index, conv_new_data_num, choose_new_data_num_limit, get_stress=False):

  '''
  gen_cp2k_task : generate cp2k tasks based on choosed struct_index

  Args :
    cp2k_dic : dictionary
      cp2k_dic includes parameters of cp2k force calculation.
    work_dir : string
      work_dir is workding directory.
    iter_id : int
      iter_id is current iteration number.
    atom_type_dic_tot : 2-d dictionary, dim = (num of lammps systems) * (num of atom types)
      example : {1:{'O':1,'H':2,'N':3},2:{'O':1,'S':2,'N':3}}
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example : {1:192,2:90}
    struct_index : dictionary (both are int)
      example : {0:{0:[2,4,6...], 1:[2,3,4...]}, 1:{0:[3,4,6...], 1:[5,6,7...]}}
                 !                !                    !
               sys_id          task_id              traj_id
    choose_new_data_num_limit : int
      choose_new_data_num_limit is the maximum number of choosed new structures.
    get_stress : bool
      get_stress is whether we need to get stress tensor
  Returns :
    box : 3-d string list
      example : [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
                 [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
    coord : 3-d string list
      example : [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
                 [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
  '''

  #copy should be done at first, because the following operators will change it!
  cp2k_param = copy.deepcopy(cp2k_dic)

  iter_dir = ''.join((work_dir + '/', 'iter_', str(iter_id)))
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
      if ( len(choosed_index) < conv_new_data_num ):
        pass
      else:
        choosed_index_num.append(len(choosed_index))
        choosed_task.append(i)

    choosed_index_num_copy = copy.deepcopy(choosed_index_num)
    if ( sum(choosed_index_num)<choose_new_data_num_limit ):
      pass
    else:
      for i in range(len(choosed_index_num)):
        choosed_index_num[i]=int(choosed_index_num_copy[i]/sum(choosed_index_num_copy)*100)

    for i in range(len(choosed_task)):
      #data file is from model0 lammps md calculation.
      choosed_index = struct_index[key][choosed_task[i]]
      choosed_index_array = np.array(choosed_index)
      np.random.shuffle(choosed_index_array)
      new_choosed_index = list(choosed_index_array[0:choosed_index_num[i]])

      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(choosed_task[i])))
      if ( new_choosed_index != [] ):
        for j in new_choosed_index:
          box_j = []
          coord_j = []
          data_file = ''.join((lammps_sys_task_dir, '/data/data_', str(j), '.lmp'))

          line_4 = linecache.getline(data_file, 4)
          line_4_split = list_dic_op.str_split(line_4, ' ')
          line_5 = linecache.getline(data_file, 5)
          line_5_split = list_dic_op.str_split(line_5, ' ')
          line_6 = linecache.getline(data_file, 6)
          line_6_split = list_dic_op.str_split(line_6, ' ')
          line_7 = linecache.getline(data_file, 7)
          line_7_split = list_dic_op.str_split(line_7, ' ')
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
          atoms_type_dic = atoms_type_dic_tot[key]
          for k in range(atoms_num):
            line_jk = linecache.getline(data_file, k+10+1)
            line_jk_split = list_dic_op.str_split(line_jk, ' ')
            atom_name = list_dic_op.get_dic_keys(atoms_type_dic, int(line_jk_split[1]))[0]
            x_str = line_jk_split[2]
            y_str = line_jk_split[3]
            z_str = line_jk_split[4].strip('\n')
            coord_j.append([atom_name,x_str,y_str,z_str])

          coord.append(coord_j)

    gen_cp2kfrc_file(cp2k_param, work_dir, iter_id, key, coord, box, get_stress)

def run_cp2kfrc(work_dir, iter_id, environ_dic, proc_num):

  '''
  run_force : perform cp2k force calculation

  Args :
    work_dir : string
      work_dir is workding directory.
    iter_id : int
      iter_id is current iteration number.
  Returns :
    none
  '''

  import subprocess

  if ( 'cp2k_exe' in environ_dic.keys() ):
    cp2k_exe = environ_dic['cp2k_exe']
  else:
    print ('These is no cp2k executable file')
    exit()

  if ( 'cp2k_env_file' in environ_dic.keys() ):
    cp2k_env = environ_dic['cp2k_env_file']
  else:
    cp2k_env = 'none'

  calc_dir = ''.join((work_dir, '/iter_', str(iter_id), '/03.cp2k_calc'))
  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(calc_dir, cmd))
  for i in range(sys_num):
    cp2k_sys_dir = ''.join((calc_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" %('task_')
    task_num = len(call.call_returns_shell(cp2k_sys_dir, cmd))

    run = '''
#! /bin/bash

source %s

mpirun -np %d %s input.inp 1> cp2k.out 2> cp2k.err
''' %(cp2k_env, proc_num, cp2k_exe)

    for j in range(task_num):
      cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(j)))
      run_file=''.join((cp2k_sys_task_dir, '/run.sh'))
      with open(run_file, 'w') as f:
        f.write(run)

      subprocess.run('chmod +x run.sh', cwd=cp2k_sys_task_dir, shell=True)
      cmd = "bash -c './run.sh'"
      subprocess.run(cmd, cwd=cp2k_sys_task_dir, shell=True)

if __name__ == '__main__':
  from collections import OrderedDict
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import cp2k_run

  work_dir = '/home/lujunbo/WORK/Deepmd/CP2K_kit/water_test'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'force_eval', 'environ']
  deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  #Test gen_cp2kfrc_file function
  coord = [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
            [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
  box = [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
         [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
#  cp2k_run.gen_cp2kfrc_file(cp2k_dic, work_dir, 1, 0, coord, box)

  #Test gen_cp2k_task function
  atoms_type_dic_tot = {0:{'O':1,'H':2}}
  atoms_num_tot = {0:3}
  struct_index = OrderedDict([(0, OrderedDict([(0, [0, 1, 2, 3, 8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 53, 54, 55, 56, 57, 62, 63, 64, 65, 66, 72, 73, 74, 75, 81, 82, 83, 90, 91, 92, 99, 100]), (1, [0, 1, 2, 3, 8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 35, 36, 37, 38, 39, 44, 45, 46, 47, 48, 53, 54, 55, 56, 57, 62, 63, 64, 65, 66, 72, 73, 74, 75, 81, 82, 83, 90, 91, 92, 99, 100]), (2, [0, 1, 2, 3, 8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 35, 36, 37, 38, 39, 45, 46, 47, 48, 53, 54, 55, 56, 57, 63, 64, 65, 66, 72, 73, 74, 75, 81, 82, 83, 84, 90, 91, 92, 93, 99, 100]), (3, [0, 1, 2, 3, 8, 9, 10, 11, 12, 17, 18, 19, 20, 21, 26, 27, 28, 29, 30, 35, 36, 37, 38, 39, 45, 46, 47, 48, 53, 54, 55, 56, 57, 63, 64, 65, 66, 72, 73, 74, 75, 81, 82, 83, 84, 90, 91, 92, 93, 99, 100])]))]) 
  conv_new_data_num = 5
  choose_new_data_num_limit = 100
  cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, 0, atoms_type_dic_tot, atoms_num_tot, \
                         struct_index, conv_new_data_num, choose_new_data_num_limit)

  #Test run_cp2kfrc function
#  cp2k_run.run_cp2kfrc(work_dir, 0, environ_dic)
