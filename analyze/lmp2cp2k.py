#! /usr/env/bin python

import os
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import read_lmp
from CP2K_kit.tools import file_tools
from CP2K_kit.analyze import check_analyze
from CP2K_kit.lib.geometry_mod import geometry

def lmp2cp2k(work_dir, lmp_log_file, lmp_traj_file, lmp_unit, atom_label, time_step, unwrap, a_vec, b_vec, c_vec):

  '''
  lmp2cp2k: transform LAMMPS trajectory into CP2K format

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit
    lmp_log_file: string
      lmp_log_file is the lammps output file.
    lmp_traj_file: string
      lmp_traj_file is the lammps trajectory file.
    lmp_unit: string
      lmp_unit is the lammps unit format. Usually, we use metal.
    atom_label: dictionary
      atom_label defines the atom name and atom id.
      Example : {1:'O', 2:'H'}
    time_step: float
      time_step is the time step of lammps md calculation.
    unwrap: bool
      unwrap is whether we need to unwrap the system.
      In general, lammps trajectory is wrapped.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example: [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example: [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
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

  atoms_num, frames_num, start_id, end_id, each = read_lmp.lmp_traj_info(lmp_traj_file, lmp_log_file)

  log_item_line_num = file_tools.grep_line_num("'Step '", lmp_log_file, work_dir)[0]

  line = linecache.getline(lmp_log_file, log_item_line_num)
  line_split = data_op.split_str(line, ' ', '\n')
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

  for i in range(frames_num):
    line = linecache.getline(lmp_log_file, i+log_item_line_num+1)
    line_split = data_op.split_str(line, ' ', '\n')
    if ( line_split[0] == 'WARNING:' ):
      log_info.log_error('There is WARNING in output file. The MD may crash, please check!')
      exit()      
    else:
      if ( 'Step' in log_item_id.keys() ):
        step.append(int(line_split[log_item_id['Step']]))
      if ( 'Temp' in log_item_id.keys() ):
        temp.append(line_split[log_item_id['Temp']])
      if ( 'PotEng' in log_item_id.keys() ):
        pot_e.append(float(line_split[log_item_id['PotEng']])*ene_lmp2cp2k)
      if ( 'KinEng' in log_item_id.keys() ):
        kin_e.append(float(line_split[log_item_id['KinEng']])*ene_lmp2cp2k)

  if ( 'Step' in log_item_id.keys() and 'Temp' in log_item_id.keys()  and \
       'PotEng' in log_item_id.keys() and 'KinEng' in log_item_id.keys() ):
    ene_file_name = ''.join((work_dir, '/cp2k-1.ener'))
    ene_file = open(ene_file_name, 'w')
    ene_file.write('#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]            Pot.[a.u.]\n')
    for i in range(frames_num):
      ene_file.write('%10d%20.6f%20.9f%20s%20.9f\n' %(step[i], step[i]*time_step, kin_e[i], temp[i], pot_e[i]))

  line = linecache.getline(lmp_traj_file, 9)
  line_split = data_op.split_str(line, ' ', '\n')
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
      line_split = data_op.split_str(line, ' ', '\n')
      atom_type.append(int(line_split[traj_item_id['type']]))
      atom_id.append(int(line_split[traj_item_id['id']]))

      if 'pos_xyz' in locals():
        pos_xyz_j = []
        pos_xyz_j.append(float(line_split[traj_item_id['x']]))
        pos_xyz_j.append(float(line_split[traj_item_id['y']]))
        pos_xyz_j.append(float(line_split[traj_item_id['z']]))

        pos_xyz.append(pos_xyz_j)

      if 'vel_xyz' in locals():
        vel_xyz_j = []
        vel_xyz_j.append(float(line_split[traj_item_id['vx']])*vel_lmp2cp2k)
        vel_xyz_j.append(float(line_split[traj_item_id['vy']])*vel_lmp2cp2k)
        vel_xyz_j.append(float(line_split[traj_item_id['vz']])*vel_lmp2cp2k)

        vel_xyz.append(vel_xyz_j)

      if 'frc_xyz' in locals():
        frc_xyz_j = []
        frc_xyz_j.append(float(line_split[traj_item_id['fx']])*frc_lmp2cp2k)
        frc_xyz_j.append(float(line_split[traj_item_id['fy']])*frc_lmp2cp2k)
        frc_xyz_j.append(float(line_split[traj_item_id['fz']])*frc_lmp2cp2k)

        frc_xyz.append(frc_xyz_j)

    atom_id_asc, asc_index = data_op.get_list_order(atom_id, 'ascend', True)
    atom_type_asc = data_op.reorder_list(atom_type, asc_index)
    if 'pos_xyz' in locals():
      pos_xyz_asc = data_op.reorder_list(pos_xyz, asc_index)
    if 'vel_xyz' in locals():
      vel_xyz_asc = data_op.reorder_list(vel_xyz, asc_index)
    if 'frc_xyz' in locals():
      frc_xyz_asc = data_op.reorder_list(frc_xyz, asc_index)

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

  linecache.clearcache()

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
  lmp2cp2k_run: kernel function to run lmp2cp2k

  Args :
    lmp2cp2k_param: dictionary
      lmp2cp2k_param contains keywords used in tranforming lammps file to cp2k file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns :
    none
  '''

  lmp2cp2k_param = check_analyze.check_lmp2cp2k_inp(lmp2cp2k_param)

  lmp_log_file = lmp2cp2k_param['lmp_log_file']
  lmp_traj_file = lmp2cp2k_param['lmp_traj_file']
  atom_label = lmp2cp2k_param['atom_label']
  time_step = lmp2cp2k_param['time_step']
  lmp_unit = lmp2cp2k_param['lmp_unit']
  unwrap = lmp2cp2k_param['unwrap']

  a_vec = lmp2cp2k_param['box']['A']
  b_vec = lmp2cp2k_param['box']['B']
  c_vec = lmp2cp2k_param['box']['C']

  print ('LMP2CP2K'.center(80, '*'), flush=True)
  lmp2cp2k(work_dir, lmp_log_file, lmp_traj_file, lmp_unit, atom_label, time_step, unwrap, a_vec, b_vec, c_vec)

  str_print = 'The CP2K format files from lammps output file and trajectory file are written in %s' %(work_dir)
  print (data_op.str_wrap(str_print, 80), flush=True)

