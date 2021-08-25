#!/usr/bin/env python

import os
import linecache
import numpy as np
from CP2K_kit.tools import atom
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import data_op
from CP2K_kit.lib import geometry_mod
from CP2K_kit.analyze import geometry
from CP2K_kit.analyze import check_analyze

def center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, center_type, center_id, \
           traj_coord_file, work_dir, file_name, trans_type=1, group_atom_1_id=[[]], group_atoms_mass=[[]]):

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
    traj_coord_file : string
      traj_coord_file is the name of trajectory file used to center.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    trans_type : int
      trans_type is the translation type when doing center.
      If trans_type is 0, we will consider the connectivity of molecules in the group.
      If trans_type is 1, we will never consider the connectivity of molecules in the group.
    group_atom_1_id : 2d int list
      group_atom_1_id is the id of first atoms in the molecules in the group.
    group_atoms_mass : 2d float list
      group_atoms_mass contains the atoms mass for each group.
  Return :
    none
  '''

  #In general, the center function will not consider connectivity.

  mass_list_len = []
  id_list_len = []
  for i in range(len(group_atoms_mass)):
    mass_list_len.append(len(group_atoms_mass[i]))
    id_list_len.append(len(group_atom_1_id[i]))

  group_atom_1_id = data_op.expand_2d_list(group_atom_1_id, max(id_list_len), 0)
  group_atoms_mass = data_op.expand_2d_list(group_atoms_mass, max(mass_list_len), 0)

  center_file_name = ''.join((work_dir, '/', file_name))
  center_file = open(center_file_name, 'w')

  for i in range(frames_num):
    #Dump atoms, atoms_mass and coordinates from trajectory file.
    coord = np.asfortranarray(np.zeros((atoms_num, 3)),dtype='float32')
    line_1 = linecache.getline(traj_coord_file, i*(atoms_num+base)+1)
    line_2 = linecache.getline(traj_coord_file, i*(atoms_num+base)+2)
    atoms = []
    atoms_mass = []
    for j in range(atoms_num):
      line_ij = linecache.getline(traj_coord_file, i*(atoms_num+base)+j+1+base+pre_base)
      line_ij_split = data_op.str_split(line_ij, ' ')
      atoms.append(line_ij_split[0])
      atoms_mass.append(atom.get_atom_mass(line_ij_split[0])[1])
      coord[j,0] = float(line_ij_split[1])
      coord[j,1] = float(line_ij_split[2])
      coord[j,2] = float(line_ij_split[3].strip('\n'))

    if ( center_type == "center_box" ):
      new_coord = geometry_mod.geometry.periodic_center_box(coord, np.asfortranarray(a_vec, dtype='float32'), \
                                                            np.asfortranarray(b_vec, dtype='float32'), \
                                                            np.asfortranarray(c_vec, dtype='float32'), \
                                                            trans_type, np.asfortranarray(group_atom_1_id, dtype='int32'), \
                                                            np.asfortranarray(group_atoms_mass, dtype='float32'), \
                                                            np.asfortranarray(id_list_len, dtype='int32'), \
                                                            np.asfortranarray(mass_list_len, dtype='int32'))
    elif ( center_type == "center_image" ):
      #We translate the atoms in the box at first, and then image
      new_coord_box = geometry_mod.geometry.periodic_center_box(coord, np.asfortranarray(a_vec, dtype='float32'), \
                                                                np.asfortranarray(b_vec, dtype='float32'), \
                                                                np.asfortranarray(c_vec, dtype='float32'), \
                                                                trans_type, np.asfortranarray(group_atom_1_id, dtype='int32'), \
                                                                np.asfortranarray(group_atoms_mass, dtype='float32'), \
                                                                np.asfortranarray(id_list_len, dtype='int32'), \
                                                                np.asfortranarray(mass_list_len, dtype='int32'))

      center_coord = np.asfortranarray(new_coord_box[center_id-1], dtype='float32')

      new_coord = geometry_mod.geometry.periodic_center_image(new_coord_box, np.asfortranarray(a_vec, dtype='float32'), \
                                                              np.asfortranarray(b_vec, dtype='float32'), \
                                                              np.asfortranarray(c_vec, dtype='float32'), center_coord, \
                                                              trans_type, np.asfortranarray(group_atom_1_id, dtype='int32'), \
                                                              np.asfortranarray(group_atoms_mass, dtype='float32'), \
                                                              np.asfortranarray(id_list_len, dtype='int32'), \
                                                              np.asfortranarray(mass_list_len, dtype='int32'))

    new_coord = geometry_mod.geometry.trans_box_center(new_coord, np.asfortranarray(atoms_mass,dtype='float32'), \
                                                       np.asfortranarray(a_vec, dtype='float32'), \
                                                       np.asfortranarray(b_vec, dtype='float32'), \
                                                       np.asfortranarray(c_vec, dtype='float32'))

    center_file.write(line_1)
    center_file.write(line_2)

    for j in range(atoms_num):
      center_file.write('%3s%21.10f%20.10f%20.10f\n' \
                        %(atoms[j], new_coord[j,0], new_coord[j,1], new_coord[j,2]))
  center_file.close()

  return center_file_name

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

  center_param = check_analyze.check_center_inp(center_param)
  center_type = center_param['center_type']
  traj_coord_file = center_param['traj_coord_file']

  if ( center_type == 'center_image' ):
    center_id = center_param['center_atom_id']
  else:
    center_id = 0

  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
  traj_info.get_traj_info(traj_coord_file, 'coord')

  log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

  a_vec = center_param['box']['A']
  b_vec = center_param['box']['B']
  c_vec = center_param['box']['C']

  print ('CENTER'.center(80,'*'), flush=True)

  if 'connect' in center_param.keys():
    atom_id = center_param['connect']['atom_id']
    group_atom = center_param['connect']['group_atom']

    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, \
    time_step, group_atom_1_id, group_atoms_mass = \
    traj_info.get_traj_info(traj_coord_file, 'coord', group_atom, atom_id, True)

    center_file = center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, center_type, center_id, \
                         traj_coord_file, work_dir, 'center.xyz', 0, group_atom_1_id, group_atoms_mass)
  else:
    center_file = center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, \
                  center_type, center_id, traj_coord_file, work_dir, 'center.xyz')

  print (data_op.str_wrap('The centered trajectory is written in %s' %(center_file), 80), flush=True)

if __name__ == '__main__':
  from CP2K_kit.tools import traj_info
  from CP2K_kit.analyze import center
  traj_coord_file = 'test.xyz'
  groups = [['Mn','F','O','O','O']]
  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, \
  time_step, group_atom_1_id, group_atoms_num = \
  traj_info.get_traj_info(file_name, groups, True)

  boxl = 18.898
  center.center(atoms_num, base, pre_base, frames_num, boxl, \
                group_atom_1_id, group_atoms_num, traj_coord_file, work_dir, 'center_box', 0)

