#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import get_cell
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import traj_tools
from CP2K_kit.analyze import rdf
from CP2K_kit.analyze import center
from CP2K_kit.analyze import check_analyze
from CP2K_kit.lib import geometry_mod
from CP2K_kit.lib import statistic_mod

def get_coord_num(atoms, coord, a_vec, b_vec, c_vec, r_cut):

  '''
  get_coord_num: get coordination number for different atom types in a trajectory file.

  Args:
    atoms: 1-d string list
      atoms is the list of atom names.
      Example: ['O', 'H', 'H', 'O', 'H', 'H']
    coord: 2-d float list, dim = (atoms_num)*3
      coord is the coordinations of atoms.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example: [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example: [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    r_cut: float
      r_cut is the cutoff value.
  Returns:
    atoms_type: 1-d string list
      atoms_type is the list of atom types.
    coord_num_avg: 1-d int list, dim = len(atoms_type)
      coord_num_avg contains the averaged coordination number for each atom type.
  '''

  atoms_type = data_op.list_replicate(atoms)
  coord_num = []

  for i in range(len(atoms_type)):
    coord_num_i = []
    for j in range(len(atoms)):
      coord_1 = []
      coord_2 = []
      if ( atoms[j] == atoms_type[i] ):
        for k in range(len(atoms)):
          if ( j != k ):
            coord_1.append(coord[j])
            coord_2.append(coord[k])

      if ( coord_1 != [] ):
        dist = geometry_mod.geometry.calculate_distance(np.asfortranarray(coord_1, dtype='float32'), \
                                                        np.asfortranarray(coord_2, dtype='float32'), \
                                                        np.asfortranarray(a_vec, dtype='float32'), \
                                                        np.asfortranarray(b_vec, dtype='float32'), \
                                                        np.asfortranarray(c_vec, dtype='float32'))
        coord_num_tmp = data_op.list_num_stat(dist, r_cut, 'less')
        coord_num_i.append(coord_num_tmp)
    coord_num.append(coord_num_i)

  coord_num_avg = []
  for i in range(len(coord_num)):
    coord_num_tmp = sum(coord_num[i])/len(coord_num[i])
    coord_num_avg.append(coord_num_tmp)

  return atoms_type, coord_num_avg

def expand_cell(atoms_num, base, pre_base, file_name, a_vec, b_vec, c_vec, a_exp, b_exp, c_exp, work_dir):

  '''
  expand_cell : expand the cell as super cell

  Args :
   atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    a_exp, b_exp, c_exp : int
      a_exp, b_exp and c_exp are the size of super cell.
      Example : 3*3*3
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  atom = []
  coord_atom = np.asfortranarray(np.zeros((atoms_num,3)),dtype='float32')
  for i in range(atoms_num):
    line_i = linecache.getline(file_name, i+base+pre_base+1)
    line_i_split = data_op.str_split(line_i, ' ')
    coord_atom[i,0] = float(line_i_split[1])
    coord_atom[i,1] = float(line_i_split[2])
    coord_atom[i,2] = float(line_i_split[3].strip('\n'))
    atom.append(line_i_split[0])

  linecache.clearcache()

  coord_atom_exp = geometry_mod.geometry.expand_cell(np.asfortranarray(coord_atom, dtype='float32'), \
                                                     np.asfortranarray(a_vec, dtype='float32'), \
                                                     np.asfortranarray(b_vec, dtype='float32'), \
                                                     np.asfortranarray(c_vec, dtype='float32'), \
                                                     a_exp, b_exp, c_exp)

  super_cell_file_name = ''.join((work_dir, '/super_cell.xyz'))
  super_cell_file = open(super_cell_file_name, 'w')
  super_atom = atom*(a_exp*b_exp*c_exp)

  super_cell.write('%d\n' %(len(super_atom)))
  for i in range(len(super_atom)):
    super_cell_file.write('%3s%21.10f%20.10f%20.10f\n' %(super_atom[i], \
    coord_atom_exp[i], coord_atom_exp[i], coord_atom_exp[i]))


def bond_length_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step, \
                     time_step, traj_coord_file, a_vec, b_vec, c_vec, atom_1_id, atom_2_id, work_dir):

  '''
  bond_length_stat: get the bond length between two atoms over different frames.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in the trajectory.
    frames_num: int
      frames_num is the number of frames in the trajectory file.
    each: int
      each is the printing frequency of of the trajectory.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    time_step: float
      time_step is the time step of md. Its unit is fs in CP2K_kit.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    atom_1_id: int
      atom_1_id is the id for atom 1.
    atom_2_id: int
      atom_2_id is the id for atom 2.
    work_dir: string
      work_dir is working directory of CP2K_kit.
  Returns:
    time: 1-d float list
      time contains time for different frames.
    distance: 1-d float array
      distance contains distance between atom 1 and atom 2 for different frames.
    distance_avg: float
      ditance_avg is the averaged distance between atom 1 and atom 2.
  '''

  a_vec = np.asfortranarray(a_vec, dtype='float32')
  b_vec = np.asfortranarray(b_vec, dtype='float32')
  c_vec = np.asfortranarray(c_vec, dtype='float32')

  center_file = center.center(atoms_num, base, pre_base, frames_num, \
                a_vec, b_vec, c_vec, 'center_box', 0, traj_coord_file, work_dir, 'center.xyz')

  frame_stat_num = int((end_step-init_step)/each+1)
  coord_atom_1 = np.asfortranarray(np.zeros((frame_stat_num,3)),dtype='float32')
  coord_atom_2 = np.asfortranarray(np.zeros((frame_stat_num,3)),dtype='float32')
  time = []

  for i in range(frame_stat_num):
    time.append(time_step*each*i)
    line_i_1 = linecache.getline(center_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+2)+atom_1_id+base)
    line_i_1_split = data_op.str_split(line_i_1, ' ')
    coord_atom_1[i,0] = float(line_i_1_split[1])
    coord_atom_1[i,1] = float(line_i_1_split[2])
    coord_atom_1[i,2] = float(line_i_1_split[3].strip('\n'))

    line_i_2 = linecache.getline(center_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+2)+atom_2_id+base)
    line_i_2_split = data_op.str_split(line_i_2, ' ')
    coord_atom_2[i,0] = float(line_i_2_split[1])
    coord_atom_2[i,1] = float(line_i_2_split[2])
    coord_atom_2[i,2] = float(line_i_2_split[3].strip('\n'))

  linecache.clearcache()

  distance = geometry_mod.geometry.calculate_distance(coord_atom_1, coord_atom_2, a_vec, b_vec, c_vec)
  distance_avg, sigma = statistic_mod.statistic.numerical_average(distance)

  cmd = 'rm %s' %(center_file)
  call.call_simple_shell(work_dir, cmd)

  return time, distance, distance_avg, sigma

def bond_angle_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, \
                    end_step, time_step, traj_coord_file, atom_1_id, atom_2_id, atom_3_id):

  #atom_2_id is the center atom for bond angle analysis.

  '''
  bond_angle_stat: get the bond angle between three atoms over different frames.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in the trajectory file.
    frames_num: int
      frames_num is the number of frames in the trajectory file.
    each: int
      each is printing frequency of md.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    atom_1_id: int
      atom_1_id is the id for atom 1.
    atom_2_id: int
      atom_2_id is the id for atom 2.
    atom_3_id: int
      atom_3_id is the id for atom 3.
  Returns:
    time: 1-d float list
      time contains time for different frames.
    angle: 1-d float array
      angle contains angle between three atoms for different frames.
    angle_avg: float
      angle_avg is the averaged angle between three atoms.
  '''

  frame_stat_num = int((end_step-init_step)/each+1)
  coord_atom_1 = np.asfortranarray(np.zeros((frame_stat_num,3)),dtype='float32')
  coord_atom_2 = np.asfortranarray(np.zeros((frame_stat_num,3)),dtype='float32')
  coord_atom_3 = np.asfortranarray(np.zeros((frame_stat_num,3)),dtype='float32')
  time = []

  for i in range(frame_stat_num):
    time.append(time_step*i*each)
    a_atom_1 = linecache.getline(traj_coord_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+2)+atom_1_id+base)
    b_atom_1 = a_atom_1.split(' ')
    c_atom_1 = []
    for j in range(len(b_atom_1)):
      if (b_atom_1[j] != ''):
        c_atom_1.append(b_atom_1[j])
    coord_atom_1[i,0] = float(c_atom_1[1])
    coord_atom_1[i,1] = float(c_atom_1[2])
    coord_atom_1[i,2] = float(c_atom_1[3].strip('\n'))

    a_atom_2 = linecache.getline(traj_coord_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+2)+atom_2_id+base)
    b_atom_2 = a_atom_2.split(' ')
    c_atom_2 = []
    for j in range(len(b_atom_2)):
      if (b_atom_2[j] != ''):
        c_atom_2.append(b_atom_2[j])
    coord_atom_2[i,0] = float(c_atom_2[1])
    coord_atom_2[i,1] = float(c_atom_2[2])
    coord_atom_2[i,2] = float(c_atom_2[3].strip('\n'))

    a_atom_3 = linecache.getline(traj_coord_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+2)+atom_3_id+base)
    b_atom_3 = a_atom_3.split(' ')
    c_atom_3 = []
    for j in range(len(b_atom_3)):
      if (b_atom_3[j] != ''):
        c_atom_3.append(b_atom_3[j])
    coord_atom_3[i,0] = float(c_atom_3[1])
    coord_atom_3[i,1] = float(c_atom_3[2])
    coord_atom_3[i,2] = float(c_atom_3[3].strip('\n'))

  linecache.clearcache()

  angle = geometry_mod.geometry.calculate_angle(coord_atom_1, coord_atom_2, coord_atom_3)
  angle_avg, sigma = statistic_mod.statistic.numerical_average(angle)

  return time, angle, angle_avg, sigma

def order_struct(atoms_num, frames_num, base, pre_base, group_atom, atom_id, traj_coord_file, a_vec, b_vec, c_vec, work_dir, file_name):

  #This function works for small molecule where there is a center atom, and other atoms are ligands.

  '''
  order_struct: reorder the system in the trajectory file.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    frames_num: int
      frames_num is the number of frames in the trajectory file.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    group_atom: 2-d string list
      group_atom is the list of name of atoms in a group.
      Example: [['Mn','F','O','O','O']]
    atom_id: 2-d int list
      atomd_id is the list of atom id for a group.
      Example: [[1,2,3,4,...,298,299,300]]
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example: [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example: [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    new_file_name: string
      new_file_name is the new file name.
    order_list: 2-d int list
      order_list is the order.
  '''

  order_list = []

  #Get the order list from the first frame
  for i in range(len(group_atom)):
    order_list_i = []
    atom_id_i = atom_id[i]
    group_atom_i = group_atom[i]
    if ( len(atom_id_i) == len(group_atom_i) ):
      for j in atom_id_i:
        order_list_i.append(j)
    else:
      if ( len(atom_id_i) != 1 and atom_id_i[0] in atom_id_i[1:len(atom_id_i)] ):
        log_info.log_error('Order structure error: complex structure is not supported!')
        exit()
      else:
        group_atom_type, group_atom_type_num = data_op.list_replicate(group_atom_i, True)

        group_coord = []
        group_atom_id = []
        for j in range(len(group_atom_type)):
          group_coord_j = []
          group_atom_id_j = []
          for k in atom_id_i:
            line_k = linecache.getline(traj_coord_file, pre_base+base+k)
            line_k_split = data_op.str_split(line_k, ' ')
            if ( line_k_split[0] == group_atom_type[j] ):
              group_coord_j.append([float(line_k_split[1]), float(line_k_split[2]), float(line_k_split[3].strip('\n'))])
              group_atom_id_j.append(k)
          group_coord.append(group_coord_j)
          group_atom_id.append(group_atom_id_j)

        for j in range(len(group_coord[0])):
          coord_1 = []
          coord_2 = []
          order_list_i.append(group_atom_id[0][j]) #The center atom
          #######################################################
          for k in range(len(group_coord)-1):
            for l in range(len(group_coord[0+k+1])):
              coord_1.append(group_coord[0][j])
              coord_2.append(group_coord[0+k+1][l])

            dist = geometry_mod.geometry.calculate_distance(np.asfortranarray(coord_1, dtype='float32'), \
                                                            np.asfortranarray(coord_2, dtype='float32'), \
                                                            np.asfortranarray(a_vec, dtype='float32'), \
                                                            np.asfortranarray(b_vec, dtype='float32'), \
                                                            np.asfortranarray(c_vec, dtype='float32'), \
                                                            len(coord_1), 3)

            dist_ascend, dist_ascend_index = data_op.list_order(dist, 'ascend', True)
            coord_1 = []
            coord_2 = []
            for l in range(group_atom_type_num[0+k+1]):
              order_list_i.append(group_atom_id[0+k+1][dist_ascend_index[l]])
          #######################################################
          #The above block is used to dump ligand atoms.
    order_list.append(order_list_i)

  new_file_name = ''.join((work_dir, '/', file_name))
  new_traj_file = open(new_file_name, 'w')
  for i in range(frames_num):
    line_i_1 = linecache.getline(traj_coord_file, (base+atoms_num)*i+1+pre_base)
    line_i_2 = linecache.getline(traj_coord_file, (base+atoms_num)*i+2+pre_base)
    new_traj_file.write(line_i_1)
    new_traj_file.write(line_i_2)
    for j in range(len(order_list)):
      for k in order_list[j]:
        line_ijk = linecache.getline(traj_coord_file, (base+atoms_num)*i+base+k+pre_base)
        new_traj_file.write(line_ijk)

  linecache.clearcache()

  return new_file_name, order_list

def first_shell(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step, atom_type_1, \
                atom_type_2, a_vec, b_vec, c_vec, traj_coord_file, dist_first_shell, dist_conv, work_dir):

  #Before you run this function, please do a rdf, then you will know the distance of first shell.

  '''
  first_shell: get the first shell number for each frame.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in the trajectory file.
    frames_num: int
      frames_num is the number of frames in the trajectory file.
    each: int
      each is printing frequency of md.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    atom_type_1: string
      atom_type_1 is the name of atom 1.
    atom_type_2: int
      atom_type_2 is the name of atom 2.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example: [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example: [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    dist_first_shell: float
      dist_first_shell is the distance of the first shell.
    dist_conv: float
      dist_conv is the converge value for the dist_first_shell.
    work_dir: string
      work_dir is working directory of CP2K_kit.
  Returns:
    first_shell: 2-d int list
      first_shell is the list of atom id in the first shell.
    dist: 2-d float list
      dist is the list of distance of the first shell.
  '''

  #distance is 3-d float list (frames_num*(number of atom_1)*(number of atom_2)).
  distance, atom_1, atom_2 = rdf.distance(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step,
                                          atom_type_1, atom_type_2, a_vec, b_vec, c_vec, traj_coord_file, work_dir)

  dim1, dim2, dim3 = np.array(distance).shape

  first_shell = []
  dist = []
  for i in range(dim1):
    first_shell_i = []
    dist_i = []
    for j in range(dim2):
      first_shell_i_j = [atom_1[j]]
      dist_i_j = []
      for k in range(dim3):
        if ( abs((distance[i][j][k]-dist_first_shell)) < dist_conv ):
          first_shell_i_j.append(atom_2[k])
          dist_i_j.append(distance[i][j][k])
      first_shell_i.append(first_shell_i_j)
      dist_i.append(dist_i_j)

    first_shell.append(first_shell_i)
    dist.append(dist_i)

  return first_shell, dist

def order_angle(center_atom_id, sur_atom_id, frame_id, each, traj_coord_file):

  '''
  This function is complicated, it is not very general. It is used in spectrum mode analysis.
  It is mainly designed to get the nearest triangle atoms in plane.

  Args:
    center_atom_id: int
      center_atom_id is the id of center atom.
    sur_atom_id: 1-d int list
      sur_atom_id is the list of atom id of surrounding atoms.
    frame_id: int
      frame_id is the id of given frame.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
  Returns:
    order : 1-d int list
      order is the order.
  '''

  atoms_num, pre_base, base, start_frame_id = traj_tools.get_block_base(traj_coord_file, 'coord')

  sur_atom_id_o = []

  coord_atom_1 = np.asfortranarray(np.zeros((1,3)),dtype='float32')
  a_atom_1 = linecache.getline(traj_coord_file, int((frame_id-start_frame_id)/each)*(atoms_num+2)+center_atom_id+base+pre_base)
  c = data_op.str_split(a_atom_1, ' ')
  coord_atom_1[0,0] = float(c[1])
  coord_atom_1[0,1] = float(c[2])
  coord_atom_1[0,2] = float(c[3].strip('\n'))

  pattern_1 = []
  pattern_2 = []
  pattern_3 = []
  for i in range(len(sur_atom_id)):
    pattern_1.append(sur_atom_id[i])
    angle_list = []
    id_list = []
    coord_atom_2 = np.asfortranarray(np.zeros((1,3)),dtype='float32')
    a_atom_2 = linecache.getline(traj_coord_file, int((frame_id-start_frame_id)/each)*(atoms_num+2)+sur_atom_id[i]+base+pre_base)
    c = data_op.str_split(a_atom_2, ' ')
    coord_atom_2[0,0] = float(c[1])
    coord_atom_2[0,1] = float(c[2])
    coord_atom_2[0,2] = float(c[3].strip('\n'))

    for j in range(len(sur_atom_id)):
      if (j != i):
        coord_atom_3 = np.asfortranarray(np.zeros((1,3)),dtype='float32')
        a_atom_3 = linecache.getline(traj_coord_file, int((frame_id-start_frame_id)/each)*(atoms_num+2)+sur_atom_id[j]+base)
        c = data_op.str_split(a_atom_3, ' ')
        coord_atom_3[0,0] = float(c[1])
        coord_atom_3[0,1] = float(c[2])
        coord_atom_3[0,2] = float(c[3].strip('\n'))
        angle = geometry_mod.geometry.calculate_angle(coord_atom_3,coord_atom_1,coord_atom_2)
        angle_list.append(angle[0])
        id_list.append(sur_atom_id[j])

    angle_list_sort, id_list_sort = zip(*sorted(zip(angle_list, id_list)))
    pattern_2.append(id_list_sort[0])
    pattern_3.append(id_list_sort[1])

  linecache.clearcache()

  order = []
  order.append(pattern_1[0])
  order.append(pattern_2[0])

  for i in range(len(pattern_1)-2):
    index_c = i+2
    a = order[index_c-2]
    b = order[index_c-1]
    for j in range(len(pattern_1)-1):
      if (pattern_1[j+1] == b and pattern_2[j+1] == a):
        order.append(pattern_3[j+1])
      if (pattern_1[j+1] == b and pattern_3[j+1] == a):
        order.append(pattern_2[j+1])

  return order

def geometry_run(geometry_param, work_dir):

  '''
  geometry_run: the kernel function of geometry module.

  Args:
    geometry_param: dictionary
      geometry_param contains keywords used in geometry functions.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  geometry_param = check_analyze.check_geometry_inp(geometry_param)

  if ( 'coord_num' in geometry_param ):
    coord_num_param = geometry_param['coord_num']

    traj_coord_file = coord_num_param['traj_coord_file']
    r_cut = coord_num_param['r_cut']
    init_step = coord_num_param['init_step']
    end_step = coord_num_param['end_step']
    a_vec = coord_num_param['box']['A']
    b_vec = coord_num_param['box']['B']
    c_vec = coord_num_param['box']['C']

    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    center_file = center.center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, \
                                'center_box', 0, traj_coord_file, work_dir, 'center.xyz')

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Coordination number of each atom type', flush=True)

    frames_num_stat = int((end_step-init_step)/each+1)

    atoms = []
    for i in range(atoms_num):
      line_i = linecache.getline(center_file, pre_base+base+i+1)
      line_i_split = data_op.str_split(line_i, ' ')
      atoms.append(line_i_split[0])

    atoms_type = data_op.list_replicate(atoms)
    coord_num_tot = [0]*len(atoms_type)

    for i in range(frames_num_stat):
      atoms = []
      coord = []
      for j in range(atoms_num):
        line_ij = linecache.getline(center_file, (atoms_num+base)*(int((init_step-start_frame_id)/each)+i)+pre_base+base+j+1)
        line_ij_split = data_op.str_split(line_ij, ' ')
        atoms.append(line_ij_split[0])
        coord.append([float(line_ij_split[1]), float(line_ij_split[2]), float(line_ij_split[3].strip('\n'))])
      atoms_type_j, coord_num_j = get_coord_num(atoms, coord, a_vec, b_vec, c_vec, r_cut)
      for j in range(len(atoms_type)):
        coord_num_tot[j] = coord_num_tot[j] + coord_num_j[j]

    linecache.clearcache()

    for i in range(len(atoms_type)):
      print ('Atom type %s is: %d' %(atoms_type[i], int(coord_num_tot[i]/frames_num_stat)), flush=True)

    cmd = 'rm %s' %(center_file)
    call.call_simple_shell(work_dir, cmd)

  elif ( 'bond_length' in geometry_param ):
    bond_length_param = geometry_param['bond_length']

    traj_coord_file = bond_length_param['traj_coord_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    atom_pair = bond_length_param['atom_pair']
    atom_1 = atom_pair[0]
    atom_2 = atom_pair[1]

    init_step = bond_length_param['init_step']
    end_step = bond_length_param['end_step']

    a_vec = geometry_param['bond_length']['box']['A']
    b_vec = geometry_param['bond_length']['box']['B']
    c_vec = geometry_param['bond_length']['box']['C']

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Analyze bond length between %d and %d' %(atom_1, atom_2), flush=True)

    time, distance, distance_avg, sigma = \
    bond_length_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, \
    end_step, time_step, traj_coord_file, a_vec, b_vec, c_vec, atom_1, atom_2, work_dir)

    dist_file = ''.join((work_dir, '/distance.csv'))
    with open(dist_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time(fs)', 'distance(Ang)'])
      for i in range(len(distance)):
        writer.writerow([time[i], distance[i]])

    str_print = 'The file containing bond length (unit: angstrom) vs time (unit: fs) is written in %s' %(dist_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    print ("The averaged bond length is %f (A) and standard error is %f (A)" %(distance_avg, sigma), flush=True)

  elif ( 'bond_angle' in geometry_param ):
    bond_angle_param = geometry_param['bond_angle']

    traj_coord_file = bond_angle_param['traj_coord_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    atom_pair = bond_angle_param['atom_pair']
    atom_1 = atom_pair[0]
    atom_2 = atom_pair[1]
    atom_3 = atom_pair[2]

    init_step = bond_angle_param['init_step']

    end_step = bond_angle_param['end_step']

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Analyze bond angle between %d and %d and %d' %(atom_1, atom_2, atom_3), flush=True)

    time, angle, angle_avg, sigma = \
    bond_angle_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step,
    end_step, time_step, traj_coord_file, atom_1, atom_2, atom_3)

    angle_file = ''.join((work_dir, '/angle.csv'))
    with open(angle_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time(fs)', 'angle(rad)'])
      for i in range(len(angle)):
        writer.writerow([time[i], angle[i]])

    str_print = 'The file containing bond angle (unit: rad) vs time (unit: fs) is written in %s' %(angle_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    print ("The averaged angle is %f (rad) and standard error is %f (rad)" %(angle_avg, sigma), flush=True)

  elif (  'first_shell' in geometry_param ):
    first_shell_param = geometry_param['first_shell']

    traj_coord_file = first_shell_param['traj_coord_file']

    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    atom_type_pair = first_shell_param['atom_type_pair']
    atom_1 = atom_type_pair[0]
    atom_2 = atom_type_pair[1]

    first_shell_dist = first_shell_param['first_shell_dist']
    dist_conv = first_shell_param['dist_conv']
    init_step = first_shell_param['init_step']
    end_step = first_shell_param['end_step']

    a_vec = first_shell_param['box']['A']
    b_vec = first_shell_param['box']['B']
    c_vec = first_shell_param['box']['C']

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Analyze first shell between %s and %s' %(atom_1, atom_2), flush=True)

    first_shell_id, dist = first_shell(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step, \
                                       atom_1, atom_2, a_vec, b_vec, c_vec, traj_coord_file, first_shell_dist, dist_conv, work_dir)

    first_shell_file = ''.join((work_dir, '/first_shell'))
    with open(first_shell_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time', 'coord_num', 'distance', 'first_shell_id'])
      for i in range(len(first_shell_id)):
        writer.writerow([time_step*i*each, len(dist[i][0]), dist[i][0], first_shell_id[i][0]])

    str_print = 'The file containing first shell information is written in %s' %(first_shell_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

  elif ( 'choose_structure' in geometry_param ):
    choose_str_param = geometry_param['choose_structure']

    traj_coord_file = choose_str_param['traj_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    init_step = choose_str_param['init_step']
    end_step = choose_str_param['end_step']
    atom_id = choose_str_param['atom_id']

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Choose structure for user defined atom id', flush=True)

    choose_str_file = traj_tools.choose_str(atoms_num, pre_base, base, each, init_step, end_step, \
                                            start_frame_id, traj_coord_file, [atom_id], work_dir, 'choose.xyz')

    str_print = 'The choosed structure file is written in %s' %(choose_str_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

  elif ( 'order_structure' in geometry_param ):
    order_str_param = geometry_param['order_structure']

    traj_coord_file = order_str_param['traj_coord_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    atom_id = order_str_param['atom_id']
    group_atom = order_str_param['group_atom']
    a_vec = order_str_param['box']['A']
    b_vec = order_str_param['box']['B']
    c_vec = order_str_param['box']['C']

    traj_choose_file = traj_tools.choose_str(atoms_num, pre_base, base, each, start_frame_id, end_frame_id, \
                       start_frame_id, traj_coord_file, atom_id, work_dir, 'traj_choose.xyz')

    atom_id_new = []
    atom_id_len = []
    for i in range(len(atom_id)):
      atom_id_len.append(len(atom_id[i]))

    for i in range(len(atom_id_len)):
      if ( i == 0 ):
        atom_id_new.append(data_op.gen_list(1, atom_id_len[i], 1))
      elif ( i == len(atom_id_len)-1 ):
        atom_id_new.append(data_op.gen_list(sum(atom_id_len)-atom_id_len[i]+1, sum(atom_id_len), 1))
      else:
        atom_id_new.append(data_op.gen_list(sum(atom_id_len[0:i])+1, sum(atom_id_len[0:i])+atom_id_len[i], 1))

    print ('GEOMETRY'.center(80, '*'), flush=True)
    print ('Order structure for the trajectory with the connectivity')
    traj_order_file, order_list = order_struct(atoms_num, frames_num, base, pre_base, group_atom, atom_id_new,\
                                  traj_choose_file, a_vec, b_vec, c_vec, work_dir, 'traj_order.xyz')
    str_print = 'The ordered structure is written in %s' %(traj_order_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    cmd = 'rm %s' %(traj_choose_file)
    call.call_simple_shell(work_dir, cmd)

if __name__ == '__main__':
  from collections import OrderedDict
  from CP2K_kit.analyze import geometry

  center_atom_id = 1
  sur_atom_id = [183, 186, 192, 195]
  frame_id = 21091
  init_step = 10000
  each = 1
  traj_coord_file = '/home/lujunbo/WORK/test/TEST_CP2K_KIT/TEST_ANALYZE/spectrum/UO22+_aimd-pos-1.xyz'
  order = order_angle(center_atom_id, sur_atom_id, frame_id, init_step, each, traj_coord_file)
  print (order)

  exit()
  atoms_num = 250
  frames_num = 801
  base = 2
  pre_base = 0
  group_tot = [OrderedDict([('atom_id', '1-250'), ('group_atom', ['Mn', 'F', 'O', 'O', 'O'])])]
#  group_tot = [{'atom_id':'1-250','group_atom':['Mn','F','O','O','O']}]
  traj_file = '/home/lujunbo/code/github/CP2K_kit/analyze/work_dir/test-pos-1.xyz'
  a_vec = [18.898,0.0,0.0]
  b_vec = [0.0,18.898,0.0]
  c_vec = [0.0,0.0,18.898]
  geometry.order_struct(atoms_num, frames_num, base, pre_base, group_tot, traj_file, a_vec, b_vec, c_vec)

