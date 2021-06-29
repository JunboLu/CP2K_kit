#! /usr/env/bin python

import os
import copy
import math
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import *

def get_box_coord(box_file, coord_file):

  '''
  get_box_coord: get box and coord from box and coord file.

  Args:
    box_file : file including cell vectors
      general box_file is like this:
      a  xx  xx  xx
      b  xx  xx  xx
      c  xx  xx  xx
    coord_file: file including atoms and coordinates
      general coord_file is like this:
      atom_a  xx  xx  xx
      atom_b  xx  xx  xx
      atom_c  xx  xx  xx
      ......
  Returns:
    tri_cell_vec : 1-d string list, dim = 6
      tri_cell_vec contains six parameters of triclinic cell.
    atoms : 1-d string list, dim = the number of atoms
      atoms contains all atom names.
    x : 1-d string list, dim = num of atoms
      x contains x component coordinate for all atoms.
    y : 1-d string list, dim = num of atoms
      y contains y component coordinate for all atoms.
    z : 1-d string list, dim = num of atoms
      z contains z component coordinate for all atoms.
  '''

  #Get box information from box file.
  cell_vec = []
  if os.path.exists(box_file):
    for i in range(3):
      cell_vec_i = []
      line_i = linecache.getline(box_file, i+1)
      line_i_split = list_dic_op.str_split(line_i, ' ')
      cell_vec_i.append(float(line_i_split[1]))
      cell_vec_i.append(float(line_i_split[2]))
      cell_vec_i.append(float(line_i_split[3].strip('\n')))
      cell_vec.append(cell_vec_i)
  else:
    print ('Can not find %s' %(box_file))
    exit()

  #Convert cell to triclinic cell. Triclinic cell just needs six parameter.
  #Triclinic cell: Lx, Ly, Lz, xy, xz, yz

  cell_vec_arr = np.array(cell_vec)
  tri_cell_a, tri_cell_b, tri_cell_c = \
  get_cell.get_triclinic_cell(cell_vec_arr[0], cell_vec_arr[1], cell_vec_arr[2])
  tri_cell_vec = []
  a_0_str = numeric.get_as_num_string(tri_cell_a[0])
  b_1_str = numeric.get_as_num_string(tri_cell_b[1])
  c_2_str = numeric.get_as_num_string(tri_cell_c[2])
  b_0_str = numeric.get_as_num_string(tri_cell_b[0])
  c_0_str = numeric.get_as_num_string(tri_cell_c[0])
  c_1_str = numeric.get_as_num_string(tri_cell_c[1])
  tri_cell_vec.append(a_0_str)
  tri_cell_vec.append(b_1_str)
  tri_cell_vec.append(c_2_str)
  tri_cell_vec.append(b_0_str)
  tri_cell_vec.append(c_0_str)
  tri_cell_vec.append(c_1_str)

  #Get atom numbers from coord file.
  if os.path.exists(coord_file):
    whole_line_num = len(open(coord_file).readlines())
    atom_num = 0
    for i in range(whole_line_num):
      line_i = linecache.getline(coord_file, i+1)
      if ( line_i != '\n' ):
        atom_num = atom_num + 1
  else:
    print ('can not find %' % (coord_file))
    exit()

  atoms = []
  x = []
  y = []
  z = []

  #Get atom coord from coord file.
  for i in range(atom_num):
    line_i = linecache.getline(coord_file, i+1)
    line_i_split = list_dic_op.str_split(line_i, ' ')
    atoms.append(line_i_split[0])
    x_float = float(line_i_split[1])
    y_float = float(line_i_split[2])
    z_float = float(line_i_split[3].strip('\n'))
    x_float_str = numeric.get_as_num_string(x_float)
    y_float_str = numeric.get_as_num_string(y_float)
    z_float_str = numeric.get_as_num_string(z_float)
    x.append(x_float_str)
    y.append(y_float_str)
    z.append(z_float_str)

  return tri_cell_vec, atoms, x, y, z

def gen_data_file(tri_cell_vec, atoms_type_index, x, y, z, task_dir, file_name):

  '''
  gen_data_file : generate lammps data file. Lammps data file contains box and coord.

  Args:
    tri_cell_vec : 1-d string list, dim = 6
      tri_cell_vec contains six parameters (Lx, Ly, Lz, x, y, z) for triclinic cell.
    atoms_type_index : 1-d int list, dim = atom numbers
      atoms_type_index contains atom type index.
    x : 1-d string list, dim = num of atoms
      x component of all atoms
    y : 1-d string list, dim = num of atoms
      y component of all atoms
    z : 1-d string list, dim = num of atoms
      z component of all atoms
    task_dir : string
      the directory to generate data file
    file_name : string
      lammps data file name
  Returns:
    atoms_type_dic_tot : dictionary, dim = the number of lammps md systems
      example : {1:{'O':1,'H':2,'N':3},2:{'O':1,'S':2,'N':3}}
    atoms_num_tot: dictionary, dim = num of lammps systems
      example : {1:192,2:90}
  '''

  data_file = open(task_dir + '/' + file_name, 'w')
  data_file.write('\n')

  atoms_type_num = len(list_dic_op.list_replicate(atoms_type_index))

  line_1 = str(len(atoms_type_index)) + ' atoms\n'
  line_2 = str(atoms_type_num) + ' atom types\n'
  data_file.write(line_1)
  data_file.write(line_2)

  line_3 = '   0.0000000000   ' + tri_cell_vec[0] + ' xlo xhi\n'
  line_4 = '   0.0000000000   ' + tri_cell_vec[1] + ' ylo yhi\n'
  line_5 = '   0.0000000000   ' + tri_cell_vec[2] + ' zlo zhi\n'
  line_6 = '   ' + tri_cell_vec[3] + ' ' + tri_cell_vec[4] + ' ' + tri_cell_vec[5] + ' xy xz yz\n'
  data_file.write(line_3)
  data_file.write(line_4)
  data_file.write(line_5)
  data_file.write(line_6)
  data_file.write('\n')
  data_file.write('Atoms # atomic\n')
  data_file.write('\n')

  for i in range(len(atoms_type_index)):
    s_1 = '     ' + str(i+1)
    s_2 = '      ' + str(atoms_type_index[i])
    s_3 = '    ' + x[i]
    s_4 = '    ' + y[i]
    if ( i != len(atoms_type_index)-1):
      s_5 = '    ' + z[i] + '\n'
    else:
      s_5 = '    ' + z[i]
    line_i = s_1 + s_2 + s_3 + s_4 + s_5
    data_file.write(line_i)

  data_file.close()

def gen_lmpmd_task(lammps_dic, work_dir, iter_id):

  '''
  gen_lmpmd_in_file : generate lammps md paramter file (.in file)

  Args:
    lammps_param : dict
      lammps_param contains parameters for lammps.
    work_dir : string
      work_dir is workding directory.
    iter_id : int
      iter_id is current iteration number.
  Returns:
    atom_type_dic_tot : dictionary
      Example : {0:{'O':1,'H':2}, 1:{'O':1,'H':2}}, The keys stands for system.
    atoms_num_tot : dictionary
      atoms_num_tot contains number of atoms for different systems.
      Example : {0:3, 1:3}
  '''

  #copy should be done at first, because the following operators will change it!
  lammps_param = copy.deepcopy(lammps_dic)

  iter_dir = ''.join((work_dir, '/iter_', str(iter_id)))
  cmd = "mkdir %s" %('02.lammps_calc')
  call.call_simple_shell(iter_dir, cmd)

  lammps_dir = ''.join((iter_dir, '/02.lammps_calc'))

  #For md simulation, lammps will use the first model.
  #The other deep potential models will run force calculation.
  train_dir_first = ''.join((work_dir, '/iter_', str(iter_id), '/01.train/0'))

  atoms_type_dic_tot = {}
  atoms_num_tot = {}

  sys_num = 0
  for key in lammps_param:
    if 'system' in key:
      sys_num = sys_num + 1

  for i in range(sys_num):
    cmd = "mkdir %s" % (''.join(('sys_', str(i))))
    call.call_simple_shell(lammps_dir, cmd)

    lammps_sys_dir = ''.join((lammps_dir, '/', 'sys_', str(i)))

    sys = 'system' + str(i)
    box_file = lammps_param[sys]['box']
    coord_file = lammps_param[sys]['coord']
    tri_cell_vec, atoms, x, y, z = get_box_coord(box_file, coord_file)

    atoms_type = list_dic_op.list_replicate(atoms)
    atoms_type_dic = {}
    for j in range(len(atoms_type)):
      atoms_type_dic[atoms_type[j]] = j+1

    atoms_type_dic_tot[i] = atoms_type_dic
    atoms_num_tot[i] = len(atoms)

    atoms_type_index = []
    for j in range(len(atoms)):
      atoms_type_index.append(atoms_type_dic[atoms[j]])

    gen_data_file(tri_cell_vec, atoms_type_index, x, y, z, lammps_sys_dir, 'data.lmp')

    #If temp and pres have different values, we will set different tasks.
    temp = lammps_param['temp']
    pres = lammps_param['pres']

    if ( isinstance(temp, str) ):
      temp = [temp]
    if ( isinstance(pres, str) ):
      pres = [pres]

    for j in range(len(temp)):
      for k in range(len(pres)):
        lammps_param_new = copy.deepcopy(lammps_param)
        lammps_param_new['temp'] = temp[j]
        lammps_param_new['pres'] = pres[k]

        cmd = "mkdir %s" % (''.join(('task_', str(j*len(pres)+k))))
        call.call_simple_shell(lammps_sys_dir, cmd)

        lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(j*len(pres)+k)))
        md_in_file = open(''.join((lammps_sys_task_dir, '/md_in.lammps')), 'w')

        for key in lammps_param:
          if ( 'system' not in key and key != 'use_metad' and key != 'plumed_file' and key != 'md_type'):
            s_1 = 'variable        '
            s_2 = key.upper() + '          '
            s_3 = 'equal ' + str(lammps_param_new[key]) + '\n'
            line_i = ''.join((s_1, s_2, s_3))
            md_in_file.write(line_i)

        md_in_file.write('\n')

        md_in_file.write('units           metal\n')
        md_in_file.write('boundary        p p p\n')
        md_in_file.write('atom_style      atomic\n')
        md_in_file.write('\n')
        md_in_file.write('neighbor        1.0 bin\n')
        md_in_file.write('neigh_modify every 1 delay 0 check yes\n')
        md_in_file.write('\n')
        md_in_file.write('box          tilt large\n')
        md_in_file.write(''.join(('read_data       ', lammps_sys_dir, '/data.lmp\n')))
        md_in_file.write('change_box   all triclinic\n')

        for l in range(len(atoms_type)):
          atom_num, atom_mass = atom.get_atom_mass(atoms_type[l])
          line_l = ''.join(('mass            ', str(l+1), ' ', str(atom_mass),'\n'))
          md_in_file.write(line_l)

        md_in_file.write('\n')
        md_in_file.write(''.join(('pair_style      deepmd ', train_dir_first, '/frozen_model.pb\n')))
        md_in_file.write('pair_coeff\n')
        md_in_file.write('\n')

        if ( 'use_metad' in lammps_param.keys() ):
          use_metad = list_dic_op.str_to_bool(lammps_param['use_metad'])
        else:
          use_metad = False

        md_in_file.write('thermo_style    custom step temp pe ke etotal press vol lx ly lz xy xz yz\n')
        md_in_file.write('thermo          ${THERMO_FREQ}\n')
        md_in_file.write('dump            1 all custom ${DUMP_FREQ} atom.dump id type x y z fx fy fz\n')
        md_in_file.write('velocity        all create ${TEMP} 835843 dist gaussian\n')

        if ( 'md_type' in lammps_param.keys() ):
          md_type = lammps_param['md_type']
        else:
          md_type = 'nvt'

        if use_metad:
          if ( 'plumed_file' in lammps_param.keys() ):
            plumed_file = lammps_param['plumed_file']
          else:
            print ('Cannot find plumed input file, please set plumed_file')
            exit()
          md_in_file.write('fix             1 all plumed plumedfile %s outfile plumed.log\n' % (plumed_file))
          if ( md_type == 'nvt' ):
            md_in_file.write('fix             2 all %s temp ${TEMP} ${TEMP} ${TAU_T}\n' % (md_type))
          elif ( md_type == 'npt' ):
            md_in_file.write('fix             2 all %s temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n' % (md_type))
        else:
          if ( md_type == 'nvt' ):
            md_in_file.write('fix             1 all %s temp ${TEMP} ${TEMP} ${TAU_T}\n' % (md_type))
          elif ( md_type == 'npt' ):
            md_in_file.write('fix             1 all %s temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n' % (md_type))
        md_in_file.write('timestep        ${TIME_STEP}\n')
        md_in_file.write('run             ${NSTEPS}\n')

        md_in_file.close()

  return atoms_type_dic_tot, atoms_num_tot

def gen_lmpfrc_file(work_dir, iter_id, atoms_num_tot, atoms_type_dic_tot):

  '''
  gen_lmpfrc_file : generate lammps parameter (.in file) for force calculations.

  Args :
    work_dir : string
      work_dir is workding directory.
    iter_id : int
      iter_id is current iteration number.
    atom_type_dic_tot : dictionary
      Example : {0:{'O':1,'H':2}, 1:{'O':1,'H':2}}, The keys stands for system.
    atoms_num_tot : dictionary
      atoms_num_tot contains number of atoms for different systems.
      Example : {0:3, 1:3}
  Returns :
    none
  '''

  #Before we run force calculation by lammps, we have run deepmd train and lammps md.
  train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  lammps_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  #We get sys_num from atoms_num_tot
  sys_num = len(atoms_num_tot)
  #We get model_num from the number of directories in train_dir.
  model_num = len(call.call_returns_shell(train_dir, 'ls'))

  for i in range(sys_num):
    lammps_sys_dir = ''.join((lammps_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lammps_sys_dir, cmd))

    for j in range(task_num):
      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(j)))
      #log.lammps is output file of lammps.
      #atom.dump is trajectory file.
      log_file = ''.join((lammps_sys_task_dir, '/log.lammps'))
      dump_file = ''.join((lammps_sys_task_dir, '/atom.dump'))

      cmd = "grep %s %s" % ('plumed', 'md_in.lammps')
      temp = call.call_returns_shell(lammps_sys_task_dir, cmd)
      if ( temp != [] and 'plumed' in temp[0] ):
        use_metad = True
      else:
        use_metad = False

      cmd = "mkdir %s" % ('data')
      call.call_simple_shell(lammps_sys_task_dir, cmd)
      model_data_dir = ''.join((lammps_sys_task_dir, '/data'))

      #Get the number of frames.
      cmd_a = "grep -n %s %s" % ("'Lx Ly Lz Xy Xz Yz'", log_file)
      a = call.call_returns_shell(lammps_sys_task_dir, cmd_a)
      cmd_b = "grep -n %s %s" % ("'Loop time'", log_file)
      b = call.call_returns_shell(lammps_sys_task_dir, cmd_b)
      a_int = int(a[0].split(':')[0])
      b_int = int(b[0].split(':')[0])
      tot_frame = b_int-a_int-1

      box_tot = []
      for k in range(tot_frame):
        box_vec = []
        line_k = linecache.getline(log_file, a_int+k+1)
        line_k_split = list_dic_op.str_split(line_k, ' ')
        for l in range(6):
          box_param_float = float(line_k_split[l+7])
          box_param_float_str = numeric.get_as_num_string(box_param_float)
          box_vec.append(box_param_float_str)
        box_tot.append(box_vec)

      #Get atom number
      atoms_num = atoms_num_tot[i]

      x_tot = []
      y_tot = []
      z_tot = []
      type_index_tot = []
      for k in range(tot_frame):
        x_k = []
        y_k = []
        z_k = []
        type_k = []
        id_k = []
        for l in range(atoms_num):
          line_kl = linecache.getline(dump_file, (atoms_num+9)*k+l+9+1)
          line_kl_split = list_dic_op.str_split(line_kl, ' ')
          id_k.append(int(line_kl_split[0]))
          type_k.append(line_kl_split[1])
          x_float = float(line_kl_split[2])
          y_float = float(line_kl_split[3])
          z_float = float(line_kl_split[4])
          x_float_str = numeric.get_as_num_string(x_float)
          y_float_str = numeric.get_as_num_string(y_float)
          z_float_str = numeric.get_as_num_string(z_float)
          x_k.append(x_float_str)
          y_k.append(y_float_str)
          z_k.append(z_float_str)

        #The atoms order in lammps trajectory is not same with that in data file.
        #So we should reorder it.
        id_k_asc, asc_index = list_dic_op.list_order(id_k, 'ascend', True)
        type_index_tot.append(list_dic_op.order_list(type_k, asc_index))
        x_tot.append(list_dic_op.order_list(x_k, asc_index))
        y_tot.append(list_dic_op.order_list(y_k, asc_index))
        z_tot.append(list_dic_op.order_list(z_k, asc_index))

      for k in range(tot_frame):
        data_file = 'data_' + str(k) + '.lmp'
        gen_data_file(box_tot[k], type_index_tot[k], x_tot[k], y_tot[k], z_tot[k], model_data_dir, data_file)

      #We run md for the first model, so force calculations are for other models.
      for k in range(model_num):
        cmd = "mkdir %s" % (''.join(('model_', str(k))))
        call.call_simple_shell(lammps_sys_task_dir, cmd)
        model_dir = ''.join((lammps_sys_task_dir, '/model_', str(k)))
        for l in range(tot_frame):
          cmd = "mkdir %s" % (''.join(('traj_', str(l))))
          call.call_simple_shell(model_dir, cmd)
          model_traj_dir = ''.join((model_dir, '/traj_', str(l)))

          if ( not use_metad and k == 0 ):
            atom_file = open(''.join((model_traj_dir, '/atom.dump')), 'w')
            for m in range(atoms_num+9):
              line_lm = linecache.getline(''.join((lammps_sys_task_dir, '/atom.dump')), l*(atoms_num+9)+m+1)
              atom_file.write(line_lm)

          else:
            frc_in_file = open(''.join((model_traj_dir, '/frc_in.lammps')), 'w')

            frc_in_file.write('units           metal\n')
            frc_in_file.write('boundary        p p p\n')
            frc_in_file.write('atom_style      atomic\n')
            frc_in_file.write('\n')
            frc_in_file.write('neighbor        1.0 bin\n')
            frc_in_file.write('\n')
            data_l_file = ''.join((model_data_dir, '/data_', str(l), '.lmp'))
            frc_in_file.write(''.join(('read_data       ', data_l_file,'\n')))

            atoms_type_i = atoms_type_dic_tot[i]
            for key in atoms_type_i:
              atom_num, atom_mass = atom.get_atom_mass(key)
              line_key = ''.join(('mass            ', str(atoms_type_i[key]), ' ', str(atom_mass),'\n'))
              frc_in_file.write(line_key)

            frc_in_file.write('\n')
            model_file = ''.join((train_dir, '/', str(k), '/frozen_model.pb'))
            frc_in_file.write(''.join(('pair_style      deepmd ', model_file, '\n')))
            frc_in_file.write('pair_coeff\n')
            frc_in_file.write('\n')

            frc_in_file.write('thermo_style    custom step temp pe ke etotal\n')
            frc_in_file.write('thermo          1\n')
            frc_in_file.write('dump            1 all custom 1 atom.dump id type x y z fx fy fz\n')
            frc_in_file.write('\n')

            frc_in_file.write('run 0\n')

            frc_in_file.close()

def run_lmpmd(work_dir, iter_id, proc_num):

  '''
  rum_lmpmd : kernel function to run lammps md.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    iter_id : int
      iter_id is the iteration id.
    proc_num : int
      proc_num is the process number of lammps.
  Returns :
    none
  '''

  import subprocess
  lammps_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(lammps_dir, cmd))

  for i in range(sys_num):
    lammps_sys_dir = ''.join((lammps_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lammps_sys_dir, cmd))

    for j in range(task_num):
      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(j)))

      cmd = "mpirun -np %d lmp < ./md_in.lammps > out" %(proc_num)
      subprocess.run(cmd, shell=True, cwd=lammps_sys_task_dir)

def lmpfrc_parallel(model_dir, parallel_num, start, end, parallel_exe):

  '''
  lmpfrc_parallel : run lammps force calculation in parallel.

  Args :
    model_dir : string
      model_dir is directory for any model.
    parallel_num : int
      parallel_num is the number of thread in parallel.
    start : int
      start is the starting number.
    end : int
      end is the endding number.
    parallel_exe : string
      parallel_exe is the parallel exacutable file.
  Returns:
    none
  '''

  import subprocess

  #run lammps in 1 thread. Here we just run force, it is a single
  #point calculation.

  run = '''
#! /bin/bash

direc=%s
parallel_num=%d
run_start=%d
run_end=%d
parallel_exe=%s

produce() {
x=$1
direc=$2
cd $direc/traj_$x
lmp < ./frc_in.lammps > out
}

export -f produce

seq $run_start $run_end | $parallel_exe -j $parallel_num produce {} $direc
''' %(model_dir, parallel_num, start, end, parallel_exe)

  run_file = ''.join((model_dir, '/run.sh'))
  with open(run_file, 'w') as f:
    f.write(run)

  subprocess.run('chmod +x run.sh', cwd=model_dir, shell=True)
  subprocess.run("bash -c './run.sh'", cwd=model_dir, shell=True)

def run_lmpfrc(work_dir, iter_id, parallel_exe, proc_num):

  '''
  rum_lmpfrc : kernel function to run lammps force calculation.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    iter_id : int
      iter_id is the iteration id.
    parallel_exe : string
      parallel_exe is parallel exacutable file.
    proc_num : int
      proc_num is the process number of lammps.
  Returns :
    none
  '''

  import subprocess
  lammps_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = 'ls | grep %s' % ('sys_')
  sys_num = len(call.call_returns_shell(lammps_dir, cmd))

  for i in range(sys_num):
    lammps_sys_dir = ''.join((lammps_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lammps_sys_dir, cmd))

    for j in range(task_num):
      lammps_sys_task_dir = ''.join((lammps_sys_dir, '/task_', str(j)))

      cmd = "grep %s %s" % ('plumed', 'md_in.lammps')
      temp = call.call_returns_shell(lammps_sys_task_dir, cmd)
      if ( temp != [] and 'plumed' in temp[0] ):
        use_metad = True
      else:
        use_metad = False

      cmd = "ls | grep %s" % ('model_')
      model_num = len(call.call_returns_shell(lammps_sys_task_dir, cmd))

      for k in range(model_num):
        model_dir = ''.join((lammps_sys_task_dir, '/model_', str(k)))

        cmd = "ls | grep %s" % ('traj_')
        traj_num = len(call.call_returns_shell(model_dir, cmd))

        start = 0
        end = start+proc_num-1
        cycle = math.ceil(traj_num/proc_num)

        parallel_num = proc_num
        for l in range(cycle):
          if ( use_metad and k == 0 ):
            lmpfrc_parallel(model_dir, parallel_num, start, end, parallel_exe)
          if ( k != 0 ):
            lmpfrc_parallel(model_dir, parallel_num, start, end, parallel_exe)
          start = start+proc_num
          end = end+proc_num
          if ( end > traj_num-1 ):
            end = traj_num-1
            parallel_num = end-start+1

if __name__ == '__main__':
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import lammps_run

  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'force_eval', 'environ']
  deepmd_dic, lammps_dic, cp2k_dic, force_eval_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)

  #Test gen_lmpmd_task function
  atoms_type_dic_tot, atoms_num_tot = \
  lammps_run.gen_lmpmd_task(lammps_dic, work_dir, 0)
  print (atoms_type_dic_tot, atoms_num_tot)

  #Test run_lmpmd function
  lammps_run.run_lmpmd(work_dir, 0, 4)

  #Test gen_lmpfrc_file function
  atoms_num_tot = {0:3}
  atoms_type_dic_tot = {0: {'O': 1, 'H': 2}}
  lammps_run.gen_lmpfrc_file(work_dir, 0, atoms_num_tot, atoms_type_dic_tot)

  #Test run_lmpfrc
  lammps_run.run_lmpfrc(work_dir, 0, 4)
