#!/usr/bin/env python

import linecache
import numpy as np
from CP2K_kit.tools import atom
from CP2K_kit.tools import traj_tools
from CP2K_kit.tools import data_op
from CP2K_kit.lib import geometry_mod

def center_of_mass(atom_name, coord):

  '''
  center_of_mass: calculate center of mass

  Args:
    atom_name: 1-d string list
      atom_name contains names of atoms in a group
      Example: ['O','H','H']
    coord: 2-d float list
      coord contains coordinate of atoms in a group
  Returns:
    center_coord: 1-d float list
  '''

  atom_num, atom_mass = atom.get_atom_mass(atom_name)

  center_coord = []

  for i in range(3):
    sum_value = 0.0
    sum_mass = 0.0
    for j in range(len(atom_name)):
      sum_value = sum_value + coord[j][i]*atom_mass[j]
      sum_mass = sum_mass + atom_mass[j]
    center_coord.append(sum_value/sum_mass)

  return center_coord

def build_droplet(work_dir, file_name, center_id, group_atom, len_conv, a_vec, b_vec, c_vec):

  '''
  build_droplet: Used to build droplet model.

  Args:
    work_dir: string
      work_dir is workding directory.
    file_name: string
      file_name is the name of supercell file.
    center_id: 1-d int list
      center_id contains centered atom id in supercell.
    group_atom: 1-d string list
      group_atom contained name of surrounded atoms.
    len_conv: float
      len_conv is the distance convergence.
    a_vec: 1-d float list
    b_vec: 1-d float list
    c_vec: 1-d float list
      a_vec, b_vec and c_vec are cell vectors.
  Returns:
    none
  '''

  coord_1 = []
  coord_2 = []

  atoms_num, pre_base, base, frame_start = traj_tools.get_block_base(file_name)

  #Get the center of mass of center part
  center_atom_coord = []
  atom_name = []
  for i in range(len(center_id)):
    line_i = linecache.getline(file_name, center_id[i]+pre_base+base)
    line_i_split = data_op.split_str(line_i, ' ', '\n')
    atom_name.append(line_i_split[0])
    center_atom_coord.append([float(line_i_split[1]), float(line_i_split[2]), float(line_i_split[3])])

  center_coord = center_of_mass(atom_name, center_atom_coord)

  all_true = []
  for i in range(len(group_atom)):
    all_true.append(True)

  index_first_atom = []
  for i in range(atoms_num-len(group_atom)):
    check_group = []
    group_atom_coord = []
    for j in range(len(group_atom)):
      line_ij = linecache.getline(file_name, pre_base+base+i+j+1)
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      group_atom_coord.append([float(line_ij_split[1]), float(line_ij_split[2]), float(line_ij_split[3])])
      if ( line_ij_split[0] == group_atom[j] ):
        check_group.append(True)
      else:
        check_group.append(False)
    if ( check_group == all_true ):
      group_center_coord = center_of_mass(group_atom, group_atom_coord)
      coord_1.append(group_center_coord)
      coord_2.append(center_coord)
      index_first_atom.append(i+1)

  coord_1_array = np.asfortranarray(coord_1, dtype='float32')
  coord_2_array = np.asfortranarray(coord_2, dtype='float32')
  a_vec_array = np.asfortranarray(a_vec, dtype='float32')
  b_vec_array = np.asfortranarray(b_vec, dtype='float32')
  c_vec_array = np.asfortranarray(c_vec, dtype='float32')
  distance = geometry_mod.geometry.calculate_distance(coord_1_array, coord_2_array, \
                                                      a_vec_array, b_vec_array, c_vec_array)

  new_file_name = ''.join((work_dir, '/', 'droplet.xyz'))
  new_file = open(new_file_name, 'w')
  for i in range(len(center_id)):
    line = linecache.getline(file_name, center_id[i]+pre_base+base)
    new_file.write(line)

  for i in range(len(index_first_atom)):
    if ( distance[i] < len_conv ):
      for j in range(len(group_atom)):
        line = linecache.getline(file_name, index_first_atom[i]+pre_base+base+j)
        new_file.write(line)

  linecache.clearcache()
  new_file.close()

def droplet_run(droplet_param, work_dir):

  '''
  droplet_run : kernel function to build droplet model

  Args:
    droplet_param: dictionary
      droplet_param contains parameter used to build droplet.
    work_dir: string
      work_dir is workding directory.
  Returns:
    none
  '''

  if ( 'traj_file' in droplet_param.keys() ):
    traj_file = droplet_param['traj_file']
  else:
    print ('Cannot find trajectory file, please set traj_file')
    exit()

  if ( 'center_id' in droplet_param.keys() ):
    center_id = droplet_param['center_id']
    center_id = data_op.get_id_list(center_id)
    center_id_num = [int(x) for x in center_id]
  else:
    print ('Cannot find id of center atom, please set center_id')
    exit()

  if ( 'group_atom' in droplet_param.keys() ):
    group_atom = droplet_param['group_atom']
  else:
    print ('Cannot find group atom, please set group_atom')
    exit()

  if ( 'len_conv' in droplet_param.keys() ):
    len_conv = float(droplet_param['len_conv'])
  else:
    len_conv = 8.0

  if ( 'box' in droplet_param.keys() ):
    A_exist = 'A' in droplet_param['box'].keys()
    B_exist = 'B' in droplet_param['box'].keys()
    C_exist = 'C' in droplet_param['box'].keys()
  else:
    print ('No box information found, please set box')
    exit()

  if ( A_exist and B_exist and C_exist):
    box_A = droplet_param['box']['A']
    box_B = droplet_param['box']['B']
    box_C = droplet_param['box']['C']
  else:
    print ('Box information wrong, please check')
    exit()

  a_vec = [float(x) for x in box_A]
  b_vec = [float(x) for x in box_B]
  c_vec = [float(x) for x in box_C]

  build_droplet(work_dir, traj_file, center_id_num, group_atom, len_conv, a_vec, b_vec, c_vec)

