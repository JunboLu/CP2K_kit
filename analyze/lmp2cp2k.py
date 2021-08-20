#! /usr/env/bin python

import os
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import list_dic_op
from CP2K_kit.lib.geometry_mod import geometry

def lmp2cp2k(work_dir, lmp_log_file, lmp_traj_file, lmp_unit, atom_label, time_step, unwrap, a_vec, b_vec, c_vec):

  '''
  lmp2cp2k : transfer LAMMPS trajectory in CP2K format

  Args :
    work_dir : string
      work_dir is the working directory of CP2K_kit
    lmp_log_file : string
      lmp_log_file is the lammps output file.
    lmp_traj_file : string
      lmp_traj_file is the lammps trajectory file.
    lmp_unit : string
      lmp_unit is the lammps unit format. Usually, we use metal.
    atom_label : dictionary
      atom_label defines the atom name and atom id.
      Example : {1:'O', 2:'H'}
    time_step : float
      time_step is the time step of lammps md calculation.
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
  Returns :
    none
  '''

  #Lammps has different units, so please be careful!
  if ( lmp_unit == 'metal' ):
    ene_lmp2cp2k = 1.0/27.2114
    pos_lmp2cp2k = 1.0
    vel_lmp2cp2k = 1.8897259886/(100000/2.4188843265857)
    frc_lmp2cp2k = 1.0/27.2114/1.8897259886

  if ( lmp_unit == 'real' ):
    ene_lmp2cp2k = 1.0/627.5094
    pos_lmp2cp2k = 1.0
    vel_lmp2cp2k = 1.8897259886/(100/2.4188843265857)
    frc_lmp2cp2k = 1.0/627.5094/1.8897259886

  line = linecache.getline(lmp_traj_file, 4)
  atoms_num = int(line.strip('\n'))
  frames_num = int(len(open(lmp_traj_file).readlines())/(atoms_num+9))

  cmd_a = "grep -n %s %s" % ("'Step '", lmp_log_file)
  a = call.call_returns_shell(work_dir, cmd_a)
  a_int = int(a[0].split(':')[0])

  line = linecache.getline(lmp_log_file, a_int)
  line_split = list_dic_op.str_split(line, ' ')
  line_split[len(line_split)-1] = line_split[len(line_split)-1].strip('\n')
  log_item_num = len(line_split)
  log_item_id = OrderedDict()

  if ( 'Step' in line_split ):
    log_item_id['Step'] = line_split.index('Step')
    step = []
  if ( 'Temp' in line_split ):
    log_item_id['Temp'] = line_split.index('Temp')
    temp = []
  if ( 'PotEng' in line_split ):
    log_item_id['PotEng'] = line_split.index('PotEng')
    pot_e = []
  if ( 'KinEng' in line_split ):
    log_item_id['KinEng'] = line_split.index('KinEng')
    kin_e = []

  frames_num_true = 0
  for i in range(frames_num):
    line = linecache.getline(lmp_log_file, i+a_int+1)
    if ( line != '' ):
      frames_num_true = frames_num_true+1
      line_split = list_dic_op.str_split(line, ' ')
      if ( 'Step' in log_item_id.keys() ):
        if ( log_item_id['Step'] != log_item_num-1 ):
          step.append(int(line_split[log_item_id['Step']]))
        else:
          step.append(int(line_split[log_item_id['Step']].strip('\n')))
      if ( 'Temp' in log_item_id.keys() ):
        if ( log_item_id['Temp'] != log_item_num-1 ):
          temp.append(line_split[log_item_id['Temp']])
        else:
          temp.append(line_split[log_item_id['Temp']].strip('\n'))
      if ( 'PotEng' in log_item_id.keys() ):
        if ( log_item_id['PotEng'] != log_item_num-1 ):
          pot_e.append(float(line_split[log_item_id['PotEng']])*ene_lmp2cp2k)
        else:
          pot_e.append(float(line_split[log_item_id['PotEng']].strip('\n'))*ene_lmp2cp2k)
      if ( 'KinEng' in log_item_id.keys() ):
        if ( log_item_id['KinEng'] != log_item_num-1 ):
          kin_e.append(float(line_split[log_item_id['KinEng']])*ene_lmp2cp2k)
        else:
          kin_e.append(float(line_split[log_item_id['KinEng']].strip('\n'))*ene_lmp2cp2k)
    else:
      break

  frames_num = frames_num_true
  if ( 'Step' in log_item_id.keys() and 'Temp' in log_item_id.keys()  and \
       'PotEng' in log_item_id.keys() and 'KinEng' in log_item_id.keys() ):
    ene_file_name = ''.join((work_dir, '/cp2k-1.ener'))
    ene_file = open(ene_file_name, 'w')
    ene_file.write('#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]\n')
    for i in range(frames_num):
      ene_file.write('%10d%20.6f%20.9f%20s%20.9f\n' %(step[i], step[i]*time_step, kin_e[i], temp[i], pot_e[i]))

  line = linecache.getline(lmp_traj_file, 9)
  line_split = list_dic_op.str_split(line, ' ')
  line_split[len(line_split)-1] = line_split[len(line_split)-1].strip('\n')
  traj_item_num = len(line_split)
  traj_item_id = OrderedDict()

  line_split_part = line_split[2:len(line_split)]

  if ( 'id' in line_split_part ):
    traj_item_id['id'] = line_split_part.index('id')
  else:
    log_info.log_error('Could not find id of atom in lammps trajectory')
    exit()

  if ( 'type' in line_split_part ):
    traj_item_id['type'] = line_split_part.index('type')
  else:
    log_info.log_error('Could not find type of atom in lammps trajectory')
    exit()

  if ( 'x' in line_split_part and 'y' in line_split_part and 'z' in line_split_part ):
    traj_item_id['x'] = line_split_part.index('x')
    traj_item_id['y'] = line_split_part.index('y')
    traj_item_id['z'] = line_split_part.index('z')
    pos_file_name = ''.join((work_dir, '/cp2k-pos-1.xyz'))
    pos_file = open(pos_file_name, 'w')
    coord = []
    atoms = []

  if ( 'vx' in line_split_part and 'vy' in line_split_part and 'vz' in line_split_part ):
    traj_item_id['vx'] = line_split_part.index('vx')
    traj_item_id['vy'] = line_split_part.index('vy')
    traj_item_id['vz'] = line_split_part.index('vz')
    vel_file_name = ''.join((work_dir, '/cp2k-vel-1.xyz'))
    vel_file = open(vel_file_name, 'w')

  if ( 'fx' in line_split_part and 'fy' in line_split_part and 'fz' in line_split_part ):
    traj_item_id['fx'] = line_split_part.index('fx')
    traj_item_id['fy'] = line_split_part.index('fy')
    traj_item_id['fz'] = line_split_part.index('fz')
    frc_file_name = ''.join((work_dir, '/cp2k-frc-1.xyz'))
    frc_file = open(frc_file_name, 'w')

  for i in range(frames_num):
    line = linecache.getline(lmp_traj_file, (atoms_num+9)*i+2)
    atom_id = []
    atom_type = []
    if ( 'x' in traj_item_id.keys() and 'y' in traj_item_id.keys() and 'z' in traj_item_id.keys() ):
      pos_xyz = []
    if ( 'vx' in traj_item_id.keys() and 'vy' in traj_item_id.keys() and 'vz' in traj_item_id.keys() ):
      vel_xyz = []
    if ( 'fx' in traj_item_id.keys() and 'fy' in traj_item_id.keys() and 'fz' in traj_item_id.keys() ):
      frc_xyz = []
    for j in range(atoms_num):
      line = linecache.getline(lmp_traj_file, (atoms_num+9)*i+9+j+1)
      line_split = list_dic_op.str_split(line, ' ')
      if ( traj_item_id['type'] != traj_item_num-1 ):
        atom_type.append(int(line_split[traj_item_id['type']]))
      else:
        atom_type.append(int(line_split[traj_item_id['type']].strip('\n')))
      if ( traj_item_id['id'] != traj_item_num-1 ):
        atom_id.append(int(line_split[traj_item_id['id']]))
      else:
        atom_type.append(int(line_split[traj_item_id['id']].strip('\n')))

      if 'pos_xyz' in locals():
        pos_xyz_j = []
        if ( traj_item_id['x'] != traj_item_num-1 ):
          pos_xyz_j.append(float(line_split[traj_item_id['x']]))
        else:
          pos_xyz_j.append(float(line_split[traj_item_id['x']].strip('\n')))
        if ( traj_item_id['y'] != traj_item_num-1 ):
          pos_xyz_j.append(float(line_split[traj_item_id['y']]))
        else:
          pos_xyz_j.append(float(line_split[traj_item_id['y']].strip('\n')))
        if ( traj_item_id['z'] != traj_item_num-1 ):
          pos_xyz_j.append(float(line_split[traj_item_id['z']]))
        else:
          pos_xyz_j.append(float(line_split[traj_item_id['z']].strip('\n')))

        pos_xyz.append(pos_xyz_j)

      if 'vel_xyz' in locals():
        vel_xyz_j = []
        if ( traj_item_id['vx'] != traj_item_num-1 ):
          vel_xyz_j.append(float(line_split[traj_item_id['vx']])*vel_lmp2cp2k)
        else:
          vel_xyz_j.append(float(line_split[traj_item_id['vx']].strip('\n'))*vel_lmp2cp2k)
        if ( traj_item_id['vy'] != traj_item_num-1 ):
          vel_xyz_j.append(float(line_split[traj_item_id['vy']])*vel_lmp2cp2k)
        else:
          vel_xyz_j.append(float(line_split[traj_item_id['vy']].strip('\n'))*vel_lmp2cp2k)
        if ( traj_item_id['vz'] != traj_item_num-1 ):
          vel_xyz_j.append(float(line_split[traj_item_id['vz']])*vel_lmp2cp2k)
        else:
          vel_xyz_j.append(float(line_split[traj_item_id['vz']].strip('\n'))*vel_lmp2cp2k)

        vel_xyz.append(vel_xyz_j)

      if 'frc_xyz' in locals():
        frc_xyz_j = []
        if ( traj_item_id['fx'] != traj_item_num-1 ):
          frc_xyz_j.append(float(line_split[traj_item_id['fx']])*frc_lmp2cp2k)
        else:
          frc_xyz_j.append(float(line_split[traj_item_id['fx']].strip('\n'))*frc_lmp2cp2k)
        if ( traj_item_id['fy'] != traj_item_num-1 ):
          frc_xyz_j.append(float(line_split[traj_item_id['fy']])*frc_lmp2cp2k)
        else:
          frc_xyz_j.append(float(line_split[traj_item_id['fy']].strip('\n'))*frc_lmp2cp2k)
        if ( traj_item_id['fz'] != traj_item_num-1 ):
          frc_xyz_j.append(float(line_split[traj_item_id['fz']])*frc_lmp2cp2k)
        else:
          frc_xyz_j.append(float(line_split[traj_item_id['fz']].strip('\n'))*frc_lmp2cp2k)
        frc_xyz.append(frc_xyz_j)

    atom_id_asc, asc_index = list_dic_op.list_order(atom_id, 'ascend', True)
    atom_type_asc = list_dic_op.order_list(atom_type, asc_index)
    if 'pos_xyz' in locals():
      pos_xyz_asc = list_dic_op.order_list(pos_xyz, asc_index)
    if 'vel_xyz' in locals():
      vel_xyz_asc = list_dic_op.order_list(vel_xyz, asc_index)
    if 'frc_xyz' in locals():
      frc_xyz_asc = list_dic_op.order_list(frc_xyz, asc_index)

    coord.append(pos_xyz_asc)
    atoms_i = []
    for j in range(atoms_num):
      atoms_i.append(atom_label[atom_type_asc[j]])
    atoms.append(atoms_i)

    if 'vel_xyz_asc' in locals():
      vel_file.write('%8d\n' %(atoms_num))
      vel_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', step[i], ', time =', step[i]*time_step, ', E =', pot_e[i]))
      for j in range(atoms_num):
        vel_file.write('%3s%21.10f%20.10f%20.10f\n' %(atom_label[atom_type_asc[j]], vel_xyz_asc[j][0], vel_xyz_asc[j][1], vel_xyz_asc[j][2]))

    if 'frc_xyz_asc' in locals():
      frc_file.write('%8d\n' %(atoms_num))
      frc_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', step[i], ', time =', step[i]*time_step, ', E =', pot_e[i]))
      for j in range(atoms_num):
        frc_file.write('%3s%21.10f%20.10f%20.10f\n' %(atom_label[atom_type_asc[j]], frc_xyz_asc[j][0], frc_xyz_asc[j][1], frc_xyz_asc[j][2]))

  #Whether we should unwrap the lammps coordinate.
  if unwrap:
    if 'coord' in locals():
      coord_array = np.asfortranarray(coord, dtype='float32')
      coord = geometry.unwrap_coord(coord_array, a_vec, b_vec, c_vec)
  if 'coord' in locals():
    for i in range(frames_num):
      pos_file.write('%8d\n' %(atoms_num))
      pos_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', step[i], ', time =', step[i]*time_step, ', E =', pot_e[i]))
      for j in range(atoms_num):
        pos_file.write('%3s%21.10f%20.10f%20.10f\n' %(atoms[i][j], coord[i][j][0], coord[i][j][1], coord[i][j][2]))

  if 'ene_file' in locals():
    ene_file.close()
  if 'pos_file' in locals():
    pos_file.close()
  if 'vel_file' in locals():
    vel_file.close()
  if 'frc_file' in locals():
    frc_file.close()

def lmp2cp2k_run(lmp2cp2k_param, work_dir):

  '''
  lmp2cp2k_run : kernel function to run lmp2cp2k

  Args :
    lmp2cp2k_param : dictionary
      lmp2cp2k_param contains keywords used in tranforming lammps file to cp2k file.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'lmp_log_file' in lmp2cp2k_param.keys() ):
    lmp_log_file = lmp2cp2k_param['lmp_log_file']
    if ( os.path.exists(lmp_log_file) ):
      pass
    else:
      log_info.log_error('%s file does not exist' %(lmp_log_file))
  else:
    log_info.log_error('No lammps output file found, please set analyze/lmp2cp2k/lmp_log_file')
    exit()

  if ( 'lmp_traj_file' in lmp2cp2k_param.keys() ):
    lmp_traj_file = lmp2cp2k_param['lmp_traj_file']
    if ( os.path.exists(lmp_traj_file) ):
      pass
    else:
      log_info.log_error('%s file does not exist' %(lmp_traj_file))
  else:
    log_info.log_error('No lammps trajectory file found, please set analyze/lmp2cp2k/lmp_traj_file')
    exit()

  if ( 'atom_label' in lmp2cp2k_param.keys() ):
    atom_label = lmp2cp2k_param['atom_label']
    atom_label_dic = OrderedDict()
    for i in range (len(atom_label)):
      lable_split = list_dic_op.str_split(atom_label[i], ':')
      atom_label_dic[int(lable_split[0])] = lable_split[1]
  else:
    log_info.log_error('No atom lable found, please set analyze/lmp2cp2k/atom_lable')
    exit()

  if ( 'time_step' in lmp2cp2k_param.keys() ):
    time_step = float(lmp2cp2k_param['time_step'])
  else:
    log_info.log_error('No time step found, please set analyze/lmp2cp2k/time_step')
    exit()

  if ( 'lmp_unit' in lmp2cp2k_param.keys() ):
    lmp_unit = lmp2cp2k_param['lmp_unit']
  else:
    log_info.log_error('No lammps unit found, please set analyze/lmp2cp2k/lmp_unit')
    exit()

  if ( 'unwrap' in lmp2cp2k_param.keys() ):
    unwrap = list_dic_op.str_to_bool(lmp2cp2k_param['unwrap'])
  else:
    unwrap = False

  if ( 'box' in lmp2cp2k_param.keys() ):
    A_exist = 'A' in lmp2cp2k_param['box'].keys()
    B_exist = 'B' in lmp2cp2k_param['box'].keys()
    C_exist = 'C' in lmp2cp2k_param['box'].keys()
  else:
    log_info.log_error('No box found, please set analyze/lmp2cp2k/box')
    exit()

  if ( A_exist and B_exist and C_exist):
    box_A = lmp2cp2k_param['box']['A']
    box_B = lmp2cp2k_param['box']['B']
    box_C = lmp2cp2k_param['box']['C']
  else:
    log_info.log_error('Box setting error, please check analzye/lmp2cp2k/box')
    exit()

  a_vec = [float(x) for x in box_A]
  b_vec = [float(x) for x in box_B]
  c_vec = [float(x) for x in box_C]

  lmp2cp2k(work_dir, lmp_log_file, lmp_traj_file, lmp_unit, atom_label_dic, time_step, unwrap, a_vec, b_vec, c_vec)
