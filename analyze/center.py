#!/usr/bin/env python

import os
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import atom
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import data_op
from CP2K_kit.lib import geometry_mod
from CP2K_kit.analyze import geometry

def center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, center_type, center_id, file_name, \
           work_dir, trans_type=1, exclude_group_id=[], group_atom_1_id=[[]], group_atoms_mass=[[]]):

  #Reference literature: Computer simulation of liquids, Oxford University press.

  '''
  center : center the system in trajectory file. This function if complicated!

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    frames_num : int
      frames_num is the number of frames in trajectory file.
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    center_type : string
      center_type is the type of center. Two choices : center_box, center_image.
    center_id : int
      center_id is the id of center atom. If we use center_box, center_id is 0.
      If we use center_image, we must set it.
    file_name : string
      file_name is the name of trajectory file used to center.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    trans_type : int
      trans_type is the translation type when doing center.
      If trans_type is 0, we will consider the connectivity of molecules in the group.
      If trans_type is 1, we will never consider the connectivity of molecules in the group.
    exclude_group_id : 1d int list
      exclude_group_id is the id of atoms that have no group.
    group_atom_1_id : 2d int list
      group_atom_1_id is the id of first atoms in the molecules in the group.
    group_atoms_mass : 2d float list
      group_atoms_mass contains the atoms mass for each group.
  Return :
    none
  '''

  #In general, the center function will not consider connectivity.
  a_vec = np.asfortranarray(a_vec, dtype='float32')
  b_vec = np.asfortranarray(b_vec, dtype='float32')
  c_vec = np.asfortranarray(c_vec, dtype='float32')
  exclude_group_id = np.asfortranarray(exclude_group_id, dtype='int32')
  group_atom_1_id = np.asfortranarray(group_atom_1_id, dtype='int32')
  group_atoms_mass = np.asfortranarray(group_atoms_mass, dtype='float32')

  if ( center_type == 'center_box' ):
    center_box_file = ''.join((work_dir, '/center_box.xyz'))
    new_file = open(center_box_file, 'w')
  elif ( center_type == 'center_image' ):
    center_image_file = ''.join((work_dir, '/center_image.xyz'))
    new_file = open(center_image_file, 'w')

  for i in range(frames_num):
    #Dump atoms, atoms_mass and coordinates from trajectory file.
    coord = np.asfortranarray(np.zeros((atoms_num, 3)),dtype='float32')
    line_1 = linecache.getline(file_name, i*(atoms_num+base)+1)
    line_2 = linecache.getline(file_name, i*(atoms_num+base)+2)
    atoms = []
    atoms_mass = []
    for j in range(atoms_num):
      line_ij = linecache.getline(file_name, i*(atoms_num+base)+j+1+base+pre_base)
      line_ij_split = data_op.str_split(line_ij, ' ')
      atoms.append(line_ij_split[0])
      atoms_mass.append(atom.get_atom_mass(line_ij_split[0])[1])
      coord[j,0] = float(line_ij_split[1])
      coord[j,1] = float(line_ij_split[2])
      coord[j,2] = float(line_ij_split[3].strip('\n'))

    if ( center_type == "center_box" ):
      new_coord = geometry_mod.geometry.periodic_center_box(coord, a_vec, b_vec, c_vec, trans_type, \
                                                            exclude_group_id, group_atom_1_id, group_atoms_mass)
    elif ( center_type == "center_image" ):
      #We translate the atoms in the box at first, and then image
      new_coord_box = geometry_mod.geometry.periodic_center_box(coord, a_vec, b_vec, c_vec, trans_type, \
                                                                exclude_group_id, group_atom_1_id, group_atoms_mass)

      center_coord = np.asfortranarray(new_coord_box[center_id-1], dtype='float32')

      new_coord = geometry_mod.geometry.periodic_center_image(new_coord_box, a_vec, b_vec, c_vec, center_coord, \
                                                              trans_type, exclude_group_id, group_atom_1_id, group_atoms_mass)
    atoms_mass = np.asfortranarray(atoms_mass,dtype='float32')
    new_coord = geometry_mod.geometry.trans_box_center(new_coord, atoms_mass, a_vec, b_vec, c_vec)

    new_file.write(line_1)
    new_file.write(line_2)

    for j in range(atoms_num):
      new_file.write('  %s        %.10f        %.10f        %.10f\n' \
                     %(atoms[j], new_coord[j,0], new_coord[j,1], new_coord[j,2]))

def center_run(center_param, work_dir):

  '''
  center_run : kernel function to run center. It will call center function.

  Args :
    center_param : dictionary
      center_param contains keywords used in center function.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'type' in center_param.keys() ):
    center_type = center_param['type']
  else:
    log_info.log_error('No center type found, please set analyze/center/type')
    exit()

  if ( 'center_atom_id' in center_param.keys() ):
    center_id = int(center_param['center_atom_id'])
  else:
    center_id = 0

  if ( 'traj_file' in center_param.keys() ):
    traj_file = center_param['traj_file']
    if ( os.path.exists(traj_file) ):
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      log_info.log_error('%s file does not exist')
      exit()
  else:
    log_info.log_error('No trajectory file found, please set analyze/center/traj_file')
    exit()

  if ( 'box' in center_param.keys() ):
    A_exist = 'A' in center_param['box'].keys()
    B_exist = 'B' in center_param['box'].keys()
    C_exist = 'C' in center_param['box'].keys()
  else:
    log_info.log_error('No box found, please set analyze/center/box')
    exit()

  if ( A_exist and B_exist and C_exist):
    box_A = center_param['box']['A']
    box_B = center_param['box']['B']
    box_C = center_param['box']['C']
  else:
    log_info.log_error('Box setting error, please check analyze/center/box')
    exit()

  a_vec = [float(x) for x in box_A]
  b_vec = [float(x) for x in box_B]
  c_vec = [float(x) for x in box_C]

  if 'connect' in center_param.keys():
    group_tot = []
    group_atom = []
    for j in center_param['connect']:
      group_atom.append(center_param['connect'][j]['group_atom'])
      group_j = OrderedDict()
      for k in center_param['connect'][j]:
        group_j[k] = center_param['connect'][j][k]
      group_tot.append(group_j)

    new_file = geometry.order_struct(atoms_num, frames_num, base, pre_base, group_tot, \
                                     traj_file, a_vec, b_vec, c_vec, work_dir)

    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, \
    time_step, exclude_group_id, group_atom_1_id, group_atoms_mass = \
    traj_info.get_traj_info(new_file, group_atom, True)

    center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, center_type, center_id, \
           new_file, work_dir, 0, exclude_group_id, group_atom_1_id, group_atoms_mass)
  else:
    center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, center_type, center_id, traj_file, work_dir)

if __name__ == '__main__':
  from CP2K_kit.tools import traj_info
  from CP2K_kit.analyze import center
  file_name = 'test.xyz'
  groups = [['Mn','F','O','O','O']]
  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, \
  time_step, exclude_group_id, group_atom_1_id, group_atoms_num = \
  traj_info.get_traj_info(file_name, groups, True)

  boxl = 18.898
  center.center(atoms_num, base, pre_base, frames_num, boxl, exclude_group_id, \
                group_atom_1_id, group_atoms_num, file_name, work_dir, 'center_box', 0)

