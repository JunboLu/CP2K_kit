#! /usr/env/bin python

import os
import math
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import *

def get_box_coord(box_file_name, coord_file_name):

  '''
  get_box_coord: get box and coord from box and coord file.

  Args:
    box_file_name: file including cell vectors
      general box_file_name is like this:
      a  xx  xx  xx
      b  xx  xx  xx
      c  xx  xx  xx
    coord_file_name: file including atoms and coordinates
      general coord_file_name is like this:
      atom_a  xx  xx  xx
      atom_b  xx  xx  xx
      atom_c  xx  xx  xx
      ......
  Returns:
    tri_cell_vec: 1-d string list, dim = 6
      tri_cell_vec contains six parameters of triclinic cell.
    atoms: 1-d string list, dim = the number of atoms
      atoms contains all atom names.
    x: 1-d string list, dim = num of atoms
      x contains x component coordinate for all atoms.
    y: 1-d string list, dim = num of atoms
      y contains y component coordinate for all atoms.
    z: 1-d string list, dim = num of atoms
      z contains z component coordinate for all atoms.
  '''

  #Get box information from box file.
  whole_line_num = len(open(box_file_name).readlines())
  cell_vec_a = []
  cell_vec_b = []
  cell_vec_c = []
  for i in range(whole_line_num):
    line_i = linecache.getline(box_file_name, i+1)
    line_i_split = data_op.split_str(line_i, ' ', '\n')
    if ( len(line_i_split) != 0 ):
      if ( len(line_i_split) == 4 and data_op.eval_str(line_i_split[0]) == 0 \
           and all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in line_i_split[1:4]) ):
        if ( line_i_split[0] == 'a' ):
          cell_vec_a.append(float(line_i_split[1]))
          cell_vec_a.append(float(line_i_split[2]))
          cell_vec_a.append(float(line_i_split[3]))
        elif ( line_i_split[0] == 'b' ):
          cell_vec_b.append(float(line_i_split[1]))
          cell_vec_b.append(float(line_i_split[2]))
          cell_vec_b.append(float(line_i_split[3]))
        elif ( line_i_split[0] == 'c' ):
          cell_vec_c.append(float(line_i_split[1]))
          cell_vec_c.append(float(line_i_split[2]))
          cell_vec_c.append(float(line_i_split[3]))
      else:
        log_info.log_error('Input error: box setting error, please check or reset deepff/lammps/system/box')
        exit()

  linecache.clearcache()

  #Convert cell to triclinic cell. Triclinic cell just needs six parameter.
  #Triclinic cell: Lx, Ly, Lz, xy, xz, yz

  if ( len(cell_vec_a) != 0 and len(cell_vec_b) != 0 and len(cell_vec_c) != 0 ):
    a, b, c, alpha, beta, gamma = get_cell.get_cell_const(np.array(cell_vec_a), np.array(cell_vec_b), np.array(cell_vec_c))
    tri_cell_a, tri_cell_b, tri_cell_c = \
    get_cell.get_triclinic_cell(a, b, c, alpha, beta, gamma)
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

  whole_line_num = len(open(coord_file_name).readlines())

  #Get atom and coord from coord file.
  atoms = []
  x = []
  y = []
  z = []
  for i in range(whole_line_num):
    line_i = linecache.getline(coord_file_name, i+1)
    line_i_split = data_op.split_str(line_i, ' ', '\n')
    if ( len(line_i_split) != 0 ):
      if ( len(line_i_split) == 4 and data_op.eval_str(line_i_split[0]) == 0 and \
           all(data_op.eval_str(x) == 1 or data_op.eval_str(x) == 2 for x in line_i_split[1:4]) ):
        atoms.append(line_i_split[0])
        x_float = float(line_i_split[1])
        y_float = float(line_i_split[2])
        z_float = float(line_i_split[3])
        x_float_str = numeric.get_as_num_string(x_float)
        y_float_str = numeric.get_as_num_string(y_float)
        z_float_str = numeric.get_as_num_string(z_float)
        x.append(x_float_str)
        y.append(y_float_str)
        z.append(z_float_str)
      else:
        log_info.log_error('Input error: coordination setting error, please check or reset deepff/lammps/system/coord')
        exit()

  linecache.clearcache()

  return tri_cell_vec, atoms, x, y, z

def get_md_sys_info(lmp_dic, tot_atoms_type_dic):

  '''
  get_md_sys_info: get the system information for lammps md.

  Args:
    lmp_dic: dictionary
      lmp_dic contains parameters for lammps.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    sys_num: int
      sys_num is the number of systems.
    atoms_type_multi_sys: 2-d dictionary, dim = (num of lammps systems) * (num of atom types)
      atoms_type_multi_sys is the atoms type for multi-systems.
      example: {0:{'O':1,'H':2,'N':3},1:{'O':1,'S':2,'N':3}}
    atoms_num_tot: dictionary
      atoms_num_tot contains number of atoms for different systems.
      Example: {0:3, 1:3}
  '''

  atoms_type_multi_sys = OrderedDict()
  atoms_num_tot = OrderedDict()
  sys_num = 0
  for key in lmp_dic:
    if 'system' in key:
      sys_num = sys_num + 1

  for i in range(sys_num):
    sys = 'system' + str(i)
    box_file_name = lmp_dic[sys]['box']
    coord_file_name = lmp_dic[sys]['coord']
    tri_cell_vec, atoms, x, y, z = get_box_coord(box_file_name, coord_file_name)
    atoms_type = data_op.list_replicate(atoms)
    atoms_type_dic = OrderedDict()
    for j in atoms_type:
      if j in tot_atoms_type_dic.keys():
        atoms_type_dic[j] = tot_atoms_type_dic[j]+1
      else:
        log_info.log_error('Input error: %s atom type in system %d is not trained, please check deepff/lammps/system' %(j, i+1))
        exit()
    atoms_type_multi_sys[i] = atoms_type_dic
    atoms_num_tot[i] = len(atoms)

  return sys_num, atoms_type_multi_sys, atoms_num_tot

def gen_data_file(tri_cell_vec, atoms_type_index, x, y, z, task_dir, file_name):

  '''
  gen_data_file: generate lammps data file. Lammps data file contains box and coord.

  Args:
    tri_cell_vec: 1-d string list, dim = 6
      tri_cell_vec contains six parameters (Lx, Ly, Lz, x, y, z) for triclinic cell.
    atoms_type_index: 1-d int list, dim = atom numbers
      atoms_type_index contains atom type index.
    x: 1-d string list, dim = num of atoms
      x component of all atoms
    y: 1-d string list, dim = num of atoms
      y component of all atoms
    z: 1-d string list, dim = num of atoms
      z component of all atoms
    task_dir: string
      the directory to generate data file
    file_name: string
      lammps data file name
  Returns:
    atoms_type_multi_sys: dictionary, dim = the number of lammps md systems
      example: {0:{'O':1,'H':2,'N':3},1:{'O':1,'S':2,'N':3}}
    atoms_num_tot: dictionary, dim = num of lammps systems
      example: {1:192,2:90}
  '''

  data_file = open(task_dir + '/' + file_name, 'w')
  data_file.write('\n')

  atoms_type_num = len(data_op.list_replicate(atoms_type_index))

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

def gen_lmpmd_task(lmp_dic, work_dir, iter_id, tot_atoms_type_dic):

  '''
  gen_lmpmd_task: generate lammps md paramter file (.in file)

  Args:
    lmp_dic: dict
      lmp_dic contains parameters for lammps.
    work_dir: string
      work_dir is workding directory.
    iter_id: int
      iter_id is current iteration number.
  Returns:
    none
  '''

  iter_dir = ''.join((work_dir, '/iter_', str(iter_id)))
  lmp_dir = ''.join((iter_dir, '/02.lammps_calc'))
  if ( not os.path.exists(lmp_dir) ):
    cmd = "mkdir %s" %('02.lammps_calc')
    call.call_simple_shell(iter_dir, cmd)

  #For md simulation, lammps will use the first model.
  #The other deep potential models will run force calculation.
  train_dir_first = ''.join((work_dir, '/iter_', str(iter_id), '/01.train/0'))

  sys_num, atoms_type_multi_sys, atoms_num_tot = get_md_sys_info(lmp_dic, tot_atoms_type_dic)

  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/', 'sys_', str(i)))
    if ( not os.path.exists(lmp_sys_dir) ):
      cmd = "mkdir %s" % (''.join(('sys_', str(i))))
      call.call_simple_shell(lmp_dir, cmd)

    sys = 'system' + str(i)
    box_file_name = lmp_dic[sys]['box']
    coord_file_name = lmp_dic[sys]['coord']
    use_metad = lmp_dic[sys]['use_metad']
    md_type = lmp_dic[sys]['md_type']

    if use_metad:
      plumed_file_name = lmp_dic[sys]['plumed_file']

    tri_cell_vec, atoms, x, y, z = get_box_coord(box_file_name, coord_file_name)
    atoms_type = data_op.list_replicate(atoms)
    atoms_type_index = []
    for j in range(len(atoms)):
      atoms_type_index.append(atoms_type_multi_sys[i][atoms[j]])

    gen_data_file(tri_cell_vec, atoms_type_index, x, y, z, lmp_sys_dir, 'data.lmp')

    #If temp and pres have different values, we will set different tasks.
    temp = lmp_dic['temp']
    pres = lmp_dic['pres']
    nsteps = int(lmp_dic['nsteps'])
    write_restart_freq = int(lmp_dic['write_restart_freq'])

    if ( isinstance(temp, float) ):
      temp = [temp]
    if ( isinstance(pres, float) ):
      pres = [pres]

    for j in range(len(temp)):
      for k in range(len(pres)):
        lmp_dic['temp'] = temp[j]
        lmp_dic['pres'] = pres[k]

        task_dir_name = ''.join(('task_', str(j*len(pres)+k)))
        lmp_sys_task_dir = ''.join((lmp_sys_dir, '/', task_dir_name))
        if ( not os.path.exists(lmp_sys_task_dir) ):
          cmd = "mkdir %s" % (task_dir_name)
          call.call_simple_shell(lmp_sys_dir, cmd)

        md_in_file = open(''.join((lmp_sys_task_dir, '/md_in.lammps')), 'w')

        cmd = "ls | grep %s" %("'tmp.restart.'")
        restart_file_name = call.call_returns_shell(lmp_sys_task_dir, cmd)
        restart_steps = []
        if ( len(restart_file_name) != 0 ):
          for i in restart_file_name:
            i_split = data_op.split_str(i, '.')
            restart_steps.append(int(i_split[len(i_split)-1]))
          restart_step = max(restart_steps)

        if ( 'restart_step' in locals() ):
          lmp_dic['nsteps'] = nsteps-restart_step
        for key in lmp_dic:
          if ( 'system' not in key and key != 'write_restart_freq' ):
            s_1 = 'variable        '
            s_2 = key.upper() + '          '
            s_3 = 'equal ' + str(lmp_dic[key]) + '\n'
            line_i = ''.join((s_1, s_2, s_3))
            md_in_file.write(line_i)

        md_in_file.write('\n')

        file_label = 0
        while True:
          dump_file_name = ''.join(('atom', str(file_label), '.dump'))
          dump_file_name_abs = ''.join((lmp_sys_task_dir, '/', dump_file_name))
          if ( os.path.exists(dump_file_name_abs) and os.path.getsize(dump_file_name_abs) != 0 ):
            file_label = file_label+1
          else:
            break

        md_in_file.write('units           metal\n')
        md_in_file.write('boundary        p p p\n')
        md_in_file.write('atom_style      atomic\n')
        md_in_file.write('timestep        ${TIME_STEP}\n')
        md_in_file.write('\n')
        md_in_file.write('neighbor        1.0 bin\n')
        md_in_file.write('neigh_modify every 1 delay 0 check yes\n')
        md_in_file.write('\n')
        md_in_file.write('box          tilt large\n')
        if ( 'restart_step' in locals() ):
          md_in_file.write('read_restart  %s\n' %(''.join(('tmp.restart.', str(restart_step)))))
        else:
          md_in_file.write(''.join(('read_data       ', lmp_sys_dir, '/data.lmp\n')))
        md_in_file.write('change_box   all triclinic\n')

        for l in range(len(atoms_type)):
          atom_num, atom_mass = atom.get_atom_mass(atoms_type[l])
          line_l = ''.join(('mass            ', str(l+1), ' ', str(atom_mass),'\n'))
          md_in_file.write(line_l)

        md_in_file.write('\n')
        md_in_file.write(''.join(('pair_style      deepmd ', train_dir_first, '/frozen_model.pb\n')))
        md_in_file.write('pair_coeff\n')
        md_in_file.write('\n')

        md_in_file.write('thermo_style    custom step temp pe ke etotal press vol lx ly lz xy xz yz\n')
        md_in_file.write('thermo          ${THERMO_FREQ}\n')
        md_in_file.write('dump            1 all custom ${DUMP_FREQ} %s id type x y z fx fy fz\n' %(dump_file_name))
        if ( 'restart_step' not in locals() ):
          md_in_file.write('velocity        all create ${TEMP} 835843 dist gaussian\n')

        if use_metad:
          md_in_file.write('fix             1 all plumed plumedfile %s outfile plumed.log\n' % (plumed_file_name))
          if ( md_type == 'nvt' ):
            md_in_file.write('fix             2 all %s temp ${TEMP} ${TEMP} ${TAU_T}\n' % (md_type))
          elif ( md_type == 'npt' ):
            md_in_file.write('fix             2 all %s temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n' % (md_type))
        else:
          if ( md_type == 'nvt' ):
            md_in_file.write('fix             1 all %s temp ${TEMP} ${TEMP} ${TAU_T}\n' % (md_type))
          elif ( md_type == 'npt' ):
            md_in_file.write('fix             1 all %s temp ${TEMP} ${TEMP} ${TAU_T} iso ${PRES} ${PRES} ${TAU_P}\n' % (md_type))
        md_in_file.write('restart  %d  tmp.restart\n' %(write_restart_freq))
        md_in_file.write('run             ${NSTEPS}\n')

        md_in_file.close()

def gen_lmpfrc_file(work_dir, iter_id, atoms_num_tot, atoms_type_multi_sys):

  '''
  gen_lmpfrc_file: generate lammps parameter (.in file) for force calculations.

  Args:
    work_dir: string
      work_dir is workding directory.
    iter_id: int
      iter_id is current iteration number.
    atom_type_dic_tot: dictionary
      Example: {0:{'O':1,'H':2}, 1:{'O':1,'H':2}}, The keys stands for system.
    atoms_num_tot: dictionary
      atoms_num_tot contains number of atoms for different systems.
      Example: {0:3, 1:3}
  Returns:
    none
  '''

  #Before we run force calculation by lammps, we have run deepmd train and lammps md.
  train_dir = ''.join((work_dir, '/iter_', str(iter_id), '/01.train'))
  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  #We get sys_num from atoms_num_tot
  sys_num = len(atoms_num_tot)
  #We get model_num from the number of directories in train_dir.
  model_num = len(call.call_returns_shell(train_dir, "ls -ll |awk '/^d/ {print $NF}'"))

  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      #log.lammps is output file of lammps.
      #atom.dump is trajectory file.
      log_file_name_abs = ''.join((lmp_sys_task_dir, '/lammps.out'))
      dump_file_name_abs = ''.join((lmp_sys_task_dir, '/atom.dump'))

      cmd = "grep %s %s" % ('plumed', 'md_in.lammps')
      plumed_info = call.call_returns_shell(lmp_sys_task_dir, cmd)
      if ( len(plumed_info) != 0 and 'plumed' in plumed_info [0] ):
        use_metad = True
      else:
        use_metad = False

      model_data_dir = ''.join((lmp_sys_task_dir, '/data'))
      if ( not os.path.exists(model_data_dir) ):
        cmd = "mkdir %s" % ('data')
        call.call_simple_shell(lmp_sys_task_dir, cmd)

      #Get the number of frames.
      line_num = file_tools.grep_line_num("'Step '", log_file_name_abs, work_dir)
      if ( line_num == 0 ):
        log_info.log_error('File error: %s file error, please check' %(log_file_name_abs))
        exit()
      else:
        a_int = line_num[0]
      line_num = file_tools.grep_line_num("'Loop time'", log_file_name_abs, work_dir)
      if ( line_num == 0 ):
        whole_line_num = len(open(log_file_name_abs, 'r').readlines())
        frames_num = whole_line_num-a_int
      else:
        frames_num = line_num[0]-a_int-1

      tot_box = []
      frames_num_true = 0
      for k in range(frames_num):
        box_vec = []
        line_k = linecache.getline(log_file_name_abs, a_int+k+1)
        line_k_split = data_op.split_str(line_k, ' ')
        if ( data_op.eval_str(line_k_split[0]) == 1 ):
          frames_num_true = frames_num_true+1
          for l in range(6):
            box_param_float = float(line_k_split[l+7])
            box_param_float_str = numeric.get_as_num_string(box_param_float)
            box_vec.append(box_param_float_str)
          tot_box.append(box_vec)

      #Get atom number
      atoms_num = atoms_num_tot[i]

      x_tot = []
      y_tot = []
      z_tot = []
      type_index_tot = []
      for k in range(frames_num_true):
        x_k = []
        y_k = []
        z_k = []
        type_k = []
        id_k = []
        for l in range(atoms_num):
          line_kl = linecache.getline(dump_file_name_abs, (atoms_num+9)*k+l+9+1)
          line_kl_split = data_op.split_str(line_kl, ' ')
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
        id_k_asc, asc_index = data_op.get_list_order(id_k, 'ascend', True)
        type_index_tot.append(data_op.reorder_list(type_k, asc_index))
        x_tot.append(data_op.reorder_list(x_k, asc_index))
        y_tot.append(data_op.reorder_list(y_k, asc_index))
        z_tot.append(data_op.reorder_list(z_k, asc_index))

      for k in range(frames_num_true):
        data_file_name = 'data_' + str(k) + '.lmp'
        gen_data_file(tot_box[k], type_index_tot[k], x_tot[k], y_tot[k], z_tot[k], model_data_dir, data_file_name)

      #We run md for the first model, so force calculations are for other models.
      for k in range(model_num):
        model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))
        if ( not os.path.exists(model_dir) ):
          cmd = "mkdir %s" % (''.join(('model_', str(k))))
          call.call_simple_shell(lmp_sys_task_dir, cmd)
        for l in range(frames_num_true):
          model_traj_dir = ''.join((model_dir, '/traj_', str(l)))
          if ( not os.path.exists(model_traj_dir) ):
            cmd = "mkdir %s" % (''.join(('traj_', str(l))))
            call.call_simple_shell(model_dir, cmd)

          if ( not use_metad and k == 0 ):
            atom_file = open(''.join((model_traj_dir, '/atom.dump')), 'w')
            for m in range(atoms_num+9):
              line_lm = linecache.getline(''.join((lmp_sys_task_dir, '/atom.dump')), l*(atoms_num+9)+m+1)
              atom_file.write(line_lm)

          else:
            frc_in_file = open(''.join((model_traj_dir, '/frc_in.lammps')), 'w')

            frc_in_file.write('units           metal\n')
            frc_in_file.write('boundary        p p p\n')
            frc_in_file.write('atom_style      atomic\n')
            frc_in_file.write('\n')
            frc_in_file.write('neighbor        1.0 bin\n')
            frc_in_file.write('\n')
            data_l_file_name_abs = ''.join((model_data_dir, '/data_', str(l), '.lmp'))
            frc_in_file.write(''.join(('read_data       ', data_l_file_name_abs,'\n')))

            atoms_type_i = atoms_type_multi_sys[i]
            for key in atoms_type_i:
              atom_num, atom_mass = atom.get_atom_mass(key)
              line_key = ''.join(('mass            ', str(atoms_type_i[key]), ' ', str(atom_mass),'\n'))
              frc_in_file.write(line_key)

            frc_in_file.write('\n')
            model_file_name_abs = ''.join((train_dir, '/', str(k), '/frozen_model.pb'))
            frc_in_file.write(''.join(('pair_style      deepmd ', model_file_name_abs, '\n')))
            frc_in_file.write('pair_coeff\n')
            frc_in_file.write('\n')

            frc_in_file.write('thermo_style    custom step temp pe ke etotal\n')
            frc_in_file.write('thermo          1\n')
            frc_in_file.write('dump            1 all custom 1 atom.dump id type x y z fx fy fz\n')
            frc_in_file.write('\n')

            frc_in_file.write('run 0\n')

            frc_in_file.close()

    linecache.clearcache()

def combine_frag_traj_file(lmp_task_dir):

  '''
  combine_frag_traj_file: combine fragement lammps trajectory files.

  Args:
    lmp_task_dir: string
      lmp_task_dir is the lammps task directory.
  Returns:
    none
  '''

  cmd = "ls | grep %s" %("'atom[0-9]'")
  atom_file_name = call.call_returns_shell(lmp_task_dir, cmd)

  if ( len(atom_file_name) == 0 ):
    log_info.log_error('lammps running error: no lammps trajectory and output file found')
    exit()
  else:
    frag_traj_num = len(atom_file_name)

  if ( frag_traj_num == 1 ):
    call.call_simple_shell(lmp_task_dir, "mv atom0.dump atom.dump")
    call.call_simple_shell(lmp_task_dir, "mv lammps0.out lammps.out")
  else:
    tot_dump_file_name_abs = ''.join((lmp_task_dir, '/atom.dump'))
    tot_log_file_name_abs = ''.join((lmp_task_dir, '/lammps.out'))
    tot_dump_file = open(tot_dump_file_name_abs, 'w')
    tot_log_file = open(tot_log_file_name_abs, 'w')
    log_file_0_name_abs = ''.join((lmp_task_dir, '/lammps0.out'))
    line_num = file_tools.grep_line_num("'Step '", log_file_0_name_abs, lmp_task_dir)
    if ( line_num == 0 ):
      log_info.log_error('File error: %s file error, lammps md task does not succeed' %(log_file_0_name_abs))
      exit()
    else:
      step_line_num = line_num[0]

    for i in range(step_line_num):
      line = linecache.getline(log_file_0_name_abs, i+1)
      tot_log_file.write(line)
    linecache.clearcache()
    for i in range(frag_traj_num):
      dump_file_name_abs = ''.join((lmp_task_dir, '/atom', str(i), '.dump'))
      log_file_name_abs = ''.join((lmp_task_dir, '/lammps', str(i), '.out'))
      atoms_num, frames_num, frames_num_fic, start_id, end_id, each = \
      read_lmp.lmp_traj_info(dump_file_name_abs, log_file_name_abs, True)
      line_num = file_tools.grep_line_num("'Step '", log_file_name_abs, lmp_task_dir)
      if ( line_num == 0 ):
        log_info.log_error('File error: %s file error, lammps md task does not succeed' %(log_file_name_abs))
        exit()
      else:
        step_line_num = line_num[0]
      for j in range(frames_num_fic):
        line = linecache.getline(log_file_name_abs, step_line_num+j+1)
        tot_log_file.write(line)
      for j in range(frames_num):
        for k in range(atoms_num+9):
          line = linecache.getline(dump_file_name_abs, (atoms_num+9)*j+k+1)
          tot_dump_file.write(line)
      if ( i == frag_traj_num-1 ):
        line_num = file_tools.grep_line_num("'Loop time'", log_file_name_abs, lmp_task_dir)
        whole_line_num = len(open(log_file_name_abs, 'r').readlines())
        if ( line_num != 0 ):
          loop_line_num = line_num[0]
          for j in range(loop_line_num, whole_line_num+1, 1):
            line = linecache.getline(log_file_name_abs, j)
            tot_log_file.write(line)
      call.call_simple_shell(lmp_task_dir, 'rm %s' %(dump_file_name_abs))
      call.call_simple_shell(lmp_task_dir, 'rm %s' %(log_file_name_abs))
      linecache.clearcache()

    tot_dump_file.close()
    tot_log_file.close()

def run_lmpmd(work_dir, iter_id, lmp_path, mpi_path, lmp_mpi_num, lmp_openmp_num, device):

  '''
  rum_lmpmd: kernel function to run lammps md.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    lmp_mpi_num: int
      lmp_mpi_num is the number of mpi ranks for lammps.
    lmp_openmp_num: int
      lmp_openmp_num is the number of openmp for lammps.
    device: 1-d int list
      device is the name of gpu devices for first node.
  Returns:
    none
  '''

  import subprocess
  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = "ls | grep %s" % ('sys_')
  sys_num = len(call.call_returns_shell(lmp_dir, cmd))

  #check generating lammps tasks
  check_lmp_md_gen = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      lmp_in_file_name_abs = ''.join((lmp_sys_task_dir, '/md_in.lammps'))
      if ( os.path.exists(lmp_in_file_name_abs) and os.path.getsize(lmp_in_file_name_abs) != 0 ):
        check_lmp_md_gen.append(0)
      else:
        check_lmp_md_gen.append(1)

  if ( len(check_lmp_md_gen) !=0 and all(i == 0 for i in check_lmp_md_gen) ):
    str_print = 'Success: generate lammps molecular dynamics tasks in %s' %(lmp_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating lammps molecular dynamics tasks error, please check iteration %d' %(iter_id))
    exit()

  #run lammps md
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))

      file_label = 0
      while True:
        log_file_name = ''.join(('lammps', str(file_label), '.out'))
        dump_file_name = ''.join(('atom', str(file_label), '.dump'))
        dump_file_name_abs = ''.join((lmp_sys_task_dir, '/', dump_file_name))
        if ( os.path.exists(dump_file_name_abs) and os.path.getsize(dump_file_name_abs) != 0 ):
          file_label = file_label+1
        else:
          break

      if ( len(device) == 0 ):
        run = '''
#! /bin/bash

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=%d

mpirun -np %d lmp < ./md_in.lammps 1> %s 2> lammps.err
''' %(lmp_path, mpi_path, lmp_openmp_num, lmp_mpi_num, log_file_name)

      else:
        device_str=data_op.comb_list_2_str(device, ',')
        run = '''
#! /bin/bash

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

export CUDA_VISIBLE_DEVICES=%s
export OMP_NUM_THREADS=%d

mpirun -np %d lmp < ./md_in.lammps 1> %s 2> lammps.err
''' %(lmp_path, mpi_path, device_str, lmp_openmp_num, lmp_mpi_num, log_file_name)

      run_file_name_abs = ''.join((lmp_sys_task_dir, '/run.sh'))
      with open(run_file_name_abs, 'w') as f:
        f.write(run)

      subprocess.run('chmod +x run.sh', cwd=lmp_sys_task_dir, shell=True)
      subprocess.run("bash -c './run.sh'", cwd=lmp_sys_task_dir, shell=True)

      combine_frag_traj_file(lmp_sys_task_dir)

  #check lammps md
  check_lmp_md_run = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      dump_file_name_abs = ''.join((lmp_sys_task_dir, '/atom.dump'))
      if ( os.path.exists(dump_file_name_abs) and os.path.getsize(dump_file_name_abs) != 0 ):
        check_lmp_md_run.append(0)
      else:
        check_lmp_md_run.append(1)

  if ( len(check_lmp_md_run) != 0 and  all(i == 0 for i in check_lmp_md_run) ):
    print ('  Success: molecular dynamics calculations for %d systems by lammps' %(sys_num))
  else:
    log_info.log_error('lammps molecular dynamics error, please check iteration %d' %(iter_id))
    exit()

def lmpfrc_parallel(model_dir, parallel_num, start, end, parallel_exe, lmp_path, mpi_path):

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

lmp_path=%s
mpi_path=%s

export PATH=$lmp_path/bin:$PATH
export PATH=$mpi_path/bin:$PATH
export LD_LIBRARY_PATH=$mpi_path/lib:$LD_LIBRARY_PATH

export OMP_NUM_THREADS=1
x=$1
direc=$2
cd $direc/traj_$x
lmp < ./frc_in.lammps 1> lammps.out 2> lammps.err
}

export -f produce

seq $run_start $run_end | $parallel_exe -j $parallel_num produce {} $direc
''' %(model_dir, parallel_num, start, end, parallel_exe, lmp_path, mpi_path)

  run_file_name_abs = ''.join((model_dir, '/run.sh'))
  with open(run_file_name_abs, 'w') as f:
    f.write(run)

  subprocess.run('chmod +x run.sh', cwd=model_dir, shell=True)
  subprocess.run("bash -c './run.sh'", cwd=model_dir, shell=True)

def run_lmpfrc(work_dir, iter_id, lmp_path, mpi_path, parallel_exe, lmp_mpi_num, atoms_num_tot):

  '''
  rum_lmpfrc: kernel function to run lammps force calculation.

  Args:
    work_dir: string
      work_dir is working directory of CP2K_kit.
    iter_id: int
      iter_id is the iteration id.
    parallel_exe: string
      parallel_exe is parallel exacutable file.
    lmp_mpi_num: int
      lmp_mpi_num is the number of mpi ranks for lammps.
  Returns :
    none
  '''

  import subprocess
  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))

  cmd = 'ls | grep %s' % ('sys_')
  sys_num = len(call.call_returns_shell(lmp_dir, cmd))

  #Check generating lammps force tasks.
  check_lmp_frc_gen = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))

      cmd = "grep %s %s" % ('plumed', 'md_in.lammps')
      plumed_info = call.call_returns_shell(lmp_sys_task_dir, cmd)
      if ( len(plumed_info) != 0 and 'plumed' in plumed_info[0] ):
        use_metad = True
      else:
        use_metad = False

      cmd = "ls | grep %s" % ('model_')
      model_num = len(call.call_returns_shell(lmp_sys_task_dir, cmd))

      for k in range(model_num):
        model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))

        cmd = "ls | grep %s" % ('traj_')
        traj_num = len(call.call_returns_shell(model_dir, cmd))

        for l in range(traj_num):
          traj_dir = ''.join((model_dir, '/traj_', str(l)))
          if ( not use_metad and k == 0 ) :
            check_file_name_abs = ''.join((traj_dir, '/atom.dump'))
          else:
            check_file_name_abs = ''.join((traj_dir, '/frc_in.lammps'))
          if ( os.path.exists(check_file_name_abs) and os.path.getsize(check_file_name_abs) != 0 ):
            check_lmp_frc_gen.append(0)
          else:
            check_lmp_frc_gen.append(1)

  if ( len(check_lmp_frc_gen) != 0 and all(i == 0 for i in check_lmp_frc_gen) ):
    str_print = 'Success: generating lammps model deviation file in %s' %(lmp_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating lammps model deviation tasks error, please check iteration %d' %(iter_id))
    exit()

  #Run lammps force.
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))

      cmd = "grep %s %s" % ('plumed', 'md_in.lammps')
      plumed_info = call.call_returns_shell(lmp_sys_task_dir, cmd)
      if ( plumed_info != [] and 'plumed' in plumed_info[0] ):
        use_metad = True
      else:
        use_metad = False

      cmd = "ls | grep %s" % ('model_')
      model_num = len(call.call_returns_shell(lmp_sys_task_dir, cmd))

      for k in range(model_num):
        model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))

        cmd = "ls | grep %s" % ('traj_')
        traj_num = len(call.call_returns_shell(model_dir, cmd))

        calculated_id = 0
        for l in range(traj_num):
          dump_file_name_abs = ''.join((model_dir, '/traj_', str(l), '/atom.dump'))
          if ( os.path.exists(dump_file_name_abs) and len(open(dump_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+9 ):
            calculated_id = l
          else:
            break

        if ( calculated_id != traj_num-1 ):
          if ( calculated_id != 0 ):
            start = calculated_id-1
          else:
            start = 0
          end = start+lmp_mpi_num-1
          if ( end > traj_num-1 ):
            end = traj_num-1
          else:
            cycle = math.ceil((traj_num-start)/lmp_mpi_num)

          for l in range(cycle):
            if ( use_metad or k != 0 ):
              lmpfrc_parallel(model_dir, end-start+1, start, end, parallel_exe, lmp_path, mpi_path)
            start = start+lmp_mpi_num
            end = end+lmp_mpi_num
            if ( end > traj_num-1 ):
              end = traj_num-1
        else:
          if ( use_metad or k != 0 ):
            lmpfrc_parallel(model_dir, 1, traj_num-1, traj_num-1, parallel_exe, lmp_path, mpi_path)

  #Check lammps force calculations.
  check_lmp_frc_run = []
  for i in range(sys_num):
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))

    cmd = "ls | grep %s" % ('task_')
    task_num = len(call.call_returns_shell(lmp_sys_dir, cmd))

    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))

      cmd = "ls | grep %s" % ('model_')
      model_num = len(call.call_returns_shell(lmp_sys_task_dir, cmd))

      for k in range(model_num):
        model_dir = ''.join((lmp_sys_task_dir, '/model_', str(k)))

        cmd = "ls | grep %s" % ('traj_')
        traj_num = len(call.call_returns_shell(model_dir, cmd))

        for l in range(traj_num):
          traj_dir = ''.join((model_dir, '/traj_', str(l)))
          dump_file_name_abs = ''.join((traj_dir, '/atom.dump'))
          if ( os.path.exists(dump_file_name_abs) and len(open(dump_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+9 ):
            check_lmp_frc_run.append(0)
          else:
            check_lmp_frc_run.append(1)

  if ( len(check_lmp_frc_run) != 0 and all(i == 0 for i in check_lmp_frc_run)):
    print ('  Success: model deviation calculations for %d systems by lammps' %(sys_num), flush=True)
  else:
    log_info.log_error('lammps model deviation calculations error, please check iteration %d' %(iter_id))
    exit()

if __name__ == '__main__':
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import lammps_run
  from CP2K_kit.deepff import check_deepff

  box_file = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md/lmp_init_data/box'
  coord_file = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md/lmp_init_data/str.inc'
  tri_cell_vec, atoms, x, y, z = get_box_coord(box_file, coord_file)
  print (atoms, x, y, z)

  exit()
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  deepmd_dic, lmp_dic, cp2k_dic, model_devi_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  proc_num = 4
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  #Test gen_lmpmd_task function
  atoms_type_multi_sys, atoms_num_tot = \
  lmp_run.gen_lmpmd_task(lmp_dic, work_dir, 0)
  #print (atoms_type_multi_sys, atoms_num_tot)

  #Test run_lmpmd function
  lmp_run.run_lmpmd(work_dir, 0, 4)

  #Test gen_lmpfrc_file function
  atoms_num_tot = {0:3}
  atoms_type_multi_sys = {0: {'O': 1, 'H': 2}}
  lmp_run.gen_lmpfrc_file(work_dir, 0, atoms_num_tot, atoms_type_multi_sys)

  #Test run_lmpfrc
  lmp_run.run_lmpfrc(work_dir, 0, 4)
