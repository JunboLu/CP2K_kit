#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from collections import OrderedDict
from CP2K_kit.tools import call
from CP2K_kit.tools import get_cell
from CP2K_kit.tools import list_dic_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import traj_tools
from CP2K_kit.analyze import rdf
from CP2K_kit.analyze import center
from CP2K_kit.lib import geometry_mod
from CP2K_kit.lib import statistic_mod

def get_coord_num(atoms, coord, a_vec, b_vec, c_vec, r_cut):

  '''
  get_coord_num : get coordination number for an atom type in a trajectory file.

  Args :
    atoms : 1-d string list
      atoms is the collection of atom names.
      Example : ['O', 'H', 'H', 'O', 'H', 'H']
    coord : 2-d float list, dim = (atoms_num)*3
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    r_cut : float
      r_cut is the cutoff value.
  Returns :
    atoms_type_coord_num : dictionary
      Example : {'O':80.0, 'H':90.0}
  '''

  atoms_type = list_dic_op.list_replicate(atoms)
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
        coord_num_temp = list_dic_op.list_num_stat(dist, r_cut, 'less')
        coord_num_i.append(coord_num_temp)
    coord_num.append(coord_num_i)

  atoms_type_coord_num = OrderedDict()
  for i in range(len(coord_num)):
    coord_num_temp = sum(coord_num[i])/len(coord_num[i])
    atoms_type_coord_num[atoms_type[i]] = coord_num_temp

  return atoms_type, atoms_type_coord_num

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
    line_i_split = list_dic_op.str_split(line_i, ' ')
    coord_atom[i,0] = float(line_i_split[1])
    coord_atom[i,1] = float(line_i_split[2])
    coord_atom[i,2] = float(line_i_split[3].strip('\n'))
    atom.append(line_i_split[0])

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


def bond_length_stat(atoms_num, base, pre_base, file_start, frames_num, each, start, end, time_step,
                     file_name, a_vec, b_vec, c_vec, atom_1_id, atom_2_id, work_dir):

  '''
  bond_length_stat : get the bond length between two atoms over different frames.

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    file_start : int
      file_start is the starting frame in trajectory file.
    frames_num : int
      frames_num is the number of frames in trajectory file.
    each : int
      each is printing frequency of md.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
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
    atom_1_id : int
      atom_1_id is the id for atom 1.
    atom_2_id : int
      atom_2_id is the id for atom 2.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    time : 1-d float list
      time contains time for different frames.
    distance : 1-d float array
      distance contains distance between atom 1 and atom 2 for different frames.
    distance_avg : float
      ditance_avg is the averaged distance between atom 1 and atom 2.
  '''

  a_vec = np.asfortranarray(a_vec, dtype='float32')
  b_vec = np.asfortranarray(b_vec, dtype='float32')
  c_vec = np.asfortranarray(c_vec, dtype='float32')

  center.center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, 'center_box', 0, file_name, work_dir)

  stat_num = int((end-start)/each+1)
  new_file_name = ''.join((work_dir, '/center_box.xyz'))
  coord_atom_1 = np.asfortranarray(np.zeros((stat_num,3)),dtype='float32')
  coord_atom_2 = np.asfortranarray(np.zeros((stat_num,3)),dtype='float32')
  time = []

  for i in range(stat_num):
    time.append(time_step*each*i)
    line_i_1 = linecache.getline(new_file_name, (int((start-file_start)/each)+i)*(atoms_num+2)+atom_1_id+base)
    line_i_1_split = list_dic_op.str_split(line_i_1, ' ')
    coord_atom_1[i,0] = float(line_i_1_split[1])
    coord_atom_1[i,1] = float(line_i_1_split[2])
    coord_atom_1[i,2] = float(line_i_1_split[3].strip('\n'))

    line_i_2 = linecache.getline(new_file_name, (int((start-file_start)/each)+i)*(atoms_num+2)+atom_2_id+base)
    line_i_2_split = list_dic_op.str_split(line_i_2, ' ')
    coord_atom_2[i,0] = float(line_i_2_split[1])
    coord_atom_2[i,1] = float(line_i_2_split[2])
    coord_atom_2[i,2] = float(line_i_2_split[3].strip('\n'))

  distance = geometry_mod.geometry.calculate_distance(coord_atom_1, coord_atom_2, a_vec, b_vec, c_vec)
  distance_avg, sigma = statistic_mod.statistic.numerical_average(distance)

  cmd = 'rm center_box.xyz'
  call.call_simple_shell(work_dir, cmd)

  return time, distance, distance_avg, sigma

def bond_angle_stat(atoms_num, base, pre_base, file_start, frames_num, each, start, end,
                    time_step, file_name, atom_1_id, atom_2_id, atom_3_id):

  #atom_2_id is the center for bond angle analysis.

  '''
  bond_angle_stat : get the bond angle between three atoms over different frames.

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    file_start : int
      file_start is the starting frame in trajectory file.
    frames_num : int
      frames_num is the number of frames in trajectory file.
    each : int
      each is printing frequency of md.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    atom_1_id : int
      atom_1_id is the id for atom 1.
    atom_2_id : int
      atom_2_id is the id for atom 2.
    atom_3_id : int
      atom_3_id is the id for atom 3.
  Returns :
    time : 1-d float list
      time contains time for different frames.
    angle : 1-d float array
      angle contains angle between three atoms for different frames.
    angle_avg : float
      angle_avg is the averaged angle between three atoms.
  '''

  stat_num = int((end-start)/each+1)
  coord_atom_1 = np.asfortranarray(np.zeros((stat_num,3)),dtype='float32')
  coord_atom_2 = np.asfortranarray(np.zeros((stat_num,3)),dtype='float32')
  coord_atom_3 = np.asfortranarray(np.zeros((stat_num,3)),dtype='float32')
  time = []

  for i in range(stat_num):
    time.append(time_step*i*each)
    a_atom_1 = linecache.getline(file_name, (int((start-file_start)/each)+i)*(atoms_num+2)+atom_1_id+base)
    b_atom_1 = a_atom_1.split(' ')
    c_atom_1 = []
    for j in range(len(b_atom_1)):
      if (b_atom_1[j] != ''):
        c_atom_1.append(b_atom_1[j])
    coord_atom_1[i,0] = float(c_atom_1[1])
    coord_atom_1[i,1] = float(c_atom_1[2])
    coord_atom_1[i,2] = float(c_atom_1[3].strip('\n'))

    a_atom_2 = linecache.getline(file_name, (int((start-file_start)/each)+i)*(atoms_num+2)+atom_2_id+base)
    b_atom_2 = a_atom_2.split(' ')
    c_atom_2 = []
    for j in range(len(b_atom_2)):
      if (b_atom_2[j] != ''):
        c_atom_2.append(b_atom_2[j])
    coord_atom_2[i,0] = float(c_atom_2[1])
    coord_atom_2[i,1] = float(c_atom_2[2])
    coord_atom_2[i,2] = float(c_atom_2[3].strip('\n'))

    a_atom_3 = linecache.getline(file_name, (int((start-file_start)/each)+i)*(atoms_num+2)+atom_3_id+base)
    b_atom_3 = a_atom_3.split(' ')
    c_atom_3 = []
    for j in range(len(b_atom_3)):
      if (b_atom_3[j] != ''):
        c_atom_3.append(b_atom_3[j])
    coord_atom_3[i,0] = float(c_atom_3[1])
    coord_atom_3[i,1] = float(c_atom_3[2])
    coord_atom_3[i,2] = float(c_atom_3[3].strip('\n'))

  angle = geometry_mod.geometry.calculate_angle(coord_atom_1, coord_atom_2, coord_atom_3)
  angle_avg, sigma = statistic_mod.statistic.numerical_average(angle)

  return time, angle, angle_avg, sigma

def order_struct(atoms_num, frames_num, base, pre_base, group_tot, traj_file, a_vec, b_vec, c_vec, work_dir):

  #This function works for small molecule where there is a center atom, and other atoms are ligands.

  '''
  order_struct : reorder the system in the trajectory file.

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    frames_num : int
      frames_num is the number of frames in trajectory file.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    group_tot : dictionary
      Example : [{'atom_id':'1-250','group_atom':['Mn','F','O','O','O']}]
    traj_file : string
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
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    new_file_name : string
      new_file_name is the new file name.
  '''

  order_list = []

  for i in range(len(group_tot)):
    atom_id = group_tot[i]['atom_id']
    atom_id_list = list_dic_op.get_id_list(atom_id)
    group_atom = group_tot[i]['group_atom']
    group_atom_type, group_atom_type_num = list_dic_op.list_replicate(group_atom, True)

    group_coord = []
    group_atom_id = []
    for j in range(len(group_atom_type)):
      group_coord_j = []
      group_atom_id_j = []
      for k in atom_id_list:
        line_k = linecache.getline(traj_file, pre_base+base+k)
        line_k_split = list_dic_op.str_split(line_k, ' ')
        if ( line_k_split[0] == group_atom_type[j] ):
          group_coord_j.append([float(line_k_split[1]), float(line_k_split[2]), float(line_k_split[3].strip('\n'))])
          group_atom_id_j.append(k)
      group_coord.append(group_coord_j)
      group_atom_id.append(group_atom_id_j)

    for j in range(len(group_coord[0])):
      coord_1 = []
      coord_2 = []
      order_list.append(group_atom_id[0][j]) #The center atom
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

        dist_ascend, dist_ascend_index = list_dic_op.list_order(dist, 'ascend', True)
        coord_1 = []
        coord_2 = []
        for l in range(group_atom_type_num[0+k+1]):
          order_list.append(group_atom_id[0+k+1][dist_ascend_index[l]])
      #######################################################
      #The above block is used to dump ligand atoms.

  exclude_list = list_dic_op.gen_list(1,atoms_num,1)
  for i in order_list:
    exclude_list.remove(i)

  order_list = exclude_list + order_list
  new_file_name = ''.join((work_dir, '/traj_order.xyz'))
  new_traj_file = open(new_file_name, 'w')
  for i in range(frames_num):
    line_i_1 = linecache.getline(traj_file, (base+atoms_num)*i+1+pre_base)
    line_i_2 = linecache.getline(traj_file, (base+atoms_num)*i+2+pre_base)
    new_traj_file.write(line_i_1)
    new_traj_file.write(line_i_2)
    for j in order_list:
      line_ij = linecache.getline(traj_file, (base+atoms_num)*i+base+j+pre_base)
      new_traj_file.write(line_ij)

  return new_file_name

def first_shell(atoms_num, base, pre_base, file_start, frames_num, each, start, end, atom_type_1, atom_type_2,
                a_vec, b_vec, c_vec, file_name, dist_first_shell, converge, work_dir):

  #Before you run this function, please do a rdf, then you will know the distance of first shell.

  '''
  first_shell : get the first shell number for each frame.

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    file_start : int
      file_start is the starting frame in trajectory file.
    frames_num : int
      frames_num is the number of frames in trajectory file.
    each : int
      each is printing frequency of md.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    atom_type_1 : string
      atom_type_1 is the name of atom 1.
    atom_type_2 : int
      atom_type_2 is the name of atom 2.
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    file_name : string
      file_name is the name of trajectory file used to analyze.
    dist_first_shell : float
      dist_first_shell is the distance of the first shell.
    converge : float
      converge is the converge value for the dist_first_shell.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    first_shell : 2-d int list
    dist : 2-d float list
  '''

  #distance is 3-d float list (frames_num*(number of atom_1)*(number of atom_2)).
  distance, atom_1, atom_2 = rdf.distance(atoms_num, base, pre_base, file_start, frames_num, each, start, end,
                                          atom_type_1, atom_type_2, a_vec, b_vec, c_vec, file_name, work_dir)

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
        if ( abs((distance[i][j][k]-dist_first_shell)) < converge ):
          first_shell_i_j.append(atom_2[k])
          dist_i_j.append(distance[i][j][k])
      first_shell_i.append(first_shell_i_j)
      dist_i.append(dist_i_j)

    first_shell.append(first_shell_i)
    dist.append(dist_i)

  return first_shell, dist

def order_angle(center_atom_id, sur_atom_id, frame_id, start, file_name):

  '''
  This function is complicated, it is not very general. It is used in spectrum mode analysis.
  It is mainly designed to get the nearest triangle atoms in plane.
  '''

  atoms_num, pre_base, base, file_start = traj_tools.get_block_base(file_name)

  sur_atom_id_o = []

  coord_atom_1 = np.asfortranarray(np.zeros((1,3)),dtype='float32')
  a_atom_1 = linecache.getline(file_name, (start-file_start+frame_id)*(atoms_num+2)+center_atom_id+base)
  c = list_dic_op.str_split(a_atom_1, ' ')
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
    a_atom_2 = linecache.getline(file_name, (start-file_start+frame_id)*(atoms_num+2)+sur_atom_id[i]+base)
    c = list_dic_op.str_split(a_atom_2, ' ')
    coord_atom_2[0,0] = float(c[1])
    coord_atom_2[0,1] = float(c[2])
    coord_atom_2[0,2] = float(c[3].strip('\n'))

    for j in range(len(sur_atom_id)):
      if (j != i):
        coord_atom_3 = np.asfortranarray(np.zeros((1,3)),dtype='float32')
        a_atom_3 = linecache.getline(file_name, (start-file_start+frame_id)*(atoms_num+2)+sur_atom_id[j]+base)
        c = list_dic_op.str_split(a_atom_3, ' ')
        coord_atom_3[0,0] = float(c[1])
        coord_atom_3[0,1] = float(c[2])
        coord_atom_3[0,2] = float(c[3].strip('\n'))
        angle = geometry_mod.geometry.calculate_angle(coord_atom_3,coord_atom_1,coord_atom_2)
        angle_list.append(angle[0])
        id_list.append(sur_atom_id[j])

    angle_list_sort, id_list_sort = zip(*sorted(zip(angle_list, id_list)))
    pattern_2.append(id_list_sort[0])
    pattern_3.append(id_list_sort[1])

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
  geometry_run : the kernel function of geometry module.

  Args :
    geometry_param : dictionary
      geometry_param contains keywords used in geometry functions.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'coord_num' in geometry_param ):
    coord_num_param = geometry_param['coord_num']

    if ( 'traj_file' in coord_num_param.keys() ):
      traj_file = coord_num_param['traj_file']
    else:
      print ('No trajectory file found, please choose the traj_file')
      exit()

    if ( 'r_cut' in coord_num_param.keys() ):
      r_cut = float(coord_num_param['r_cut'])
    else:
      r_cut = 6.0

    if ( 'box' in coord_num_param.keys() ):
      A_exist = 'A' in coord_num_param['box'].keys()
      B_exist = 'B' in coord_num_param['box'].keys()
      C_exist = 'C' in coord_num_param['box'].keys()
    else:
      print ('No box information found, please set box')
      exit()

    if ( A_exist and B_exist and C_exist):
      box_A = coord_num_param['box']['A']
      box_B = coord_num_param['box']['B']
      box_C = coord_num_param['box']['C']
    else:
      print ('Box information wrong, please check')
      exit()

    a_vec = [float(x) for x in box_A]
    b_vec = [float(x) for x in box_B]
    c_vec = [float(x) for x in box_C]

    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_file)

    coord_num_file = ''.join((work_dir, '/coord_num.csv'))
    with open(coord_num_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['frame','coord_num'])
      for i in range(frames_num):
        atoms = []
        coord = []
        for j in range(atoms_num):
          line_ij = linecache.getline(traj_file, (atoms_num+base)*i+pre_base+base+j+1)
          line_ij_split = list_dic_op.str_split(line_ij, ' ')
          atoms.append(line_ij_split[0])
          coord.append([float(line_ij_split[1]), float(line_ij_split[2]), float(line_ij_split[3].strip('\n'))])
        atoms_type, coord_num = get_coord_num(atoms, coord, a_vec, b_vec, c_vec, r_cut)
        writer.writerow([i, coord_num])

  elif ( 'bond_length' in geometry_param ):
    bond_length_param = geometry_param['bond_length']

    if ( 'traj_file' in bond_length_param.keys() ):
      traj_file = bond_length_param['traj_file']
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      print ('No trajectory file found, please choose the traj_file')
      exit()

    if ( 'atom_pair' in bond_length_param.keys() ):
      atom_pair = [int (x) for x in bond_length_param['atom_pair']]
      atom_1 = atom_pair[0]
      atom_2 = atom_pair[1]
    else:
      print ('No atom pair used to calculate bond length found, please set atom_pair')
      exit()

    if ( 'init_step' in bond_length_param.keys() ):
      init_step = int(bond_length_param['init_step'])
    else:
      init_step = start_frame_id

    if ( 'end_step' in bond_length_param.keys() ):
      end_step = int(bond_length_param['end_step'])
    else:
      end_step = start_frame_id + each

    if ( 'box' in bond_length_param.keys() ):
      A_exist = 'A' in bond_length_param['box'].keys()
      B_exist = 'B' in bond_length_param['box'].keys()
      C_exist = 'C' in bond_length_param['box'].keys()
    else:
      print ('No box information found, please set box')
      exit()

    if ( A_exist and B_exist and C_exist):
      box_A = bond_length_param['box']['A']
      box_B = bond_length_param['box']['B']
      box_C = bond_length_param['box']['C']
    else:
      print ('Box information wrong, please check')
      exit()

    a_vec = [float(x) for x in box_A]
    b_vec = [float(x) for x in box_B]
    c_vec = [float(x) for x in box_C]

    time, distance, distance_avg, sigma = \
    bond_length_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, \
    end_step, time_step, traj_file, a_vec, b_vec, c_vec, atom_1, atom_2, work_dir)

    dist_file = ''.join((work_dir, '/distance.csv'))
    with open(dist_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time', 'distance'])
      for i in range(len(distance)):
        writer.writerow([time[i], distance[i]])

    print ("The averaged bond length is %f (A) and standard error is %f (A)" %(distance_avg, sigma))

  elif ( 'bond_angle' in geometry_param ):
    bond_angle_param = geometry_param['bond_angle']

    if ( 'traj_file' in bond_angle_param.keys() ):
      traj_file = bond_angle_param['traj_file']
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      print ('No trajectory file found, please choose the traj_file')
      exit()

    if ( 'atom_pair' in bond_angle_param.keys() ):
      atom_pair = [int (x) for x in bond_angle_param['atom_pair']]
      atom_1 = atom_pair[0]
      atom_2 = atom_pair[1]
      atom_3 = atom_pair[2]
    else:
      print ('No atom pair used to calculate bond length found, please set atom_pair')
      exit()

    if ( 'init_step' in bond_angle_param.keys() ):
      init_step = int(bond_angle_param['init_step'])
    else:
      init_step = start_frame_id

    if ( 'end_step' in bond_angle_param.keys() ):
      end_step = int(bond_angle_param['end_step'])
    else:
      end_step = start_frame_id + each

    time, angle, angle_avg, sigma = \
    bond_angle_stat(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step,
    end_step, time_step, traj_file, atom_1, atom_2, atom_3)

    angle_file = ''.join((work_dir, '/angle.csv'))
    with open(angle_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time', 'angle'])
      for i in range(len(angle)):
        writer.writerow([time[i], angle[i]])

    print ("The averaged angle is %f (rad) and standard error is %f (rad)" %(angle_avg, sigma))

  elif (  'first_shell' in geometry_param ):
    first_shell_param = geometry_param['first_shell']

    if ( 'traj_file' in first_shell_param.keys() ):
      traj_file = first_shell_param['traj_file']
      if ( os.path.exists(traj_file) ):
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_file)
      else:
        print ('Cannot find %s file' % (traj_file))
        exit()
    else:
      print ('No trajectory file found, please choose the traj_file')
      exit()

    if ( 'atom_type_pair' in first_shell_param.keys() ):
      atom_type_pair = first_shell_param['atom_type_pair']
      atom_1 = atom_type_pair[0]
      atom_2 = atom_type_pair[1]
    else:
      print ('No atom type pair found, please set atom_type_pair')
      exit()

    if ( 'first_shell_dist' in first_shell_param.keys() ):
      first_shell_dist = float(first_shell_param['first_shell_dist'])
    else:
      print ('No first shell distance found, please run rdf and then set first_shell_dist')
      exit()

    if ( 'dist_conv' in first_shell_param.keys() ):
      dist_conv = float(first_shell_param['dist_conv'])
    else:
      dist_conv = 0.1

    if ( 'init_step' in first_shell_param.keys() ):
      init_step = int(first_shell_param['init_step'])
    else:
      init_step = start_frame_id

    if ( 'end_step' in first_shell_param.keys() ):
      end_step = int(first_shell_param['end_step'])
    else:
      end_step = start_frame_id + each

    if ( 'box' in first_shell_param.keys() ):
      A_exist = 'A' in first_shell_param['box'].keys()
      B_exist = 'B' in first_shell_param['box'].keys()
      C_exist = 'C' in first_shell_param['box'].keys()
    else:
      print ('No box information found, please set box')
      exit()

    if ( A_exist and B_exist and C_exist):
      box_A = first_shell_param['box']['A']
      box_B = first_shell_param['box']['B']
      box_C = first_shell_param['box']['C']
    else:
      print ('Box information wrong, please check')
      exit()

    a_vec = [float(x) for x in box_A]
    b_vec = [float(x) for x in box_B]
    c_vec = [float(x) for x in box_C]

    first_shell_id, dist = first_shell(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step, \
                                       atom_1, atom_2, a_vec, b_vec, c_vec, traj_file, first_shell_dist, dist_conv, work_dir)

    first_shell_file = open('first_shell','w')
    for i in range(len(first_shell_id)):
      first_shell_file.write(str(first_shell_id[i])+' '+str(dist[i])+'\n')

  elif (  'choose_structure' in geometry_param ):
    choose_str_param = geometry_param['choose_structure']

    if ( 'traj_file' in choose_str_param.keys() ):
      traj_file = choose_str_param['traj_file']
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      print ('No trajectory file found, please choose the traj_file')
      exit()

    if ( 'init_step' in choose_str_param.keys() ):
      init_step = int(choose_str_param['init_step'])
    else:
      init_step = start_frame_id

    if ( 'end_step' in choose_str_param.keys() ):
      end_step = int(choose_str_param['end_step'])
    else:
      end_step = start_frame_id + each

    if ( 'atom_id' in choose_str_param.keys() ):
      atom_id = list_dic_op.get_id_list(choose_str_param['atom_id'])
    else:
      print ('No atom id found, please set atom_id')

    traj_tools.choose_str(atoms_num, pre_base, base, each, init_step, end_step, \
                          start_frame_id, traj_file, atom_id, work_dir)

if __name__ == '__main__':
  from collections import OrderedDict
  from CP2K_kit.analyze import geometry
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

