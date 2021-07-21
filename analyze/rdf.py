#!/usr/bin/env python

import os
import csv
import math
import linecache
import numpy as np
from CP2K_kit.tools import list_dic_op
from CP2K_kit.analyze import center
from CP2K_kit.tools import traj_info
from CP2K_kit.lib import geometry_mod

def distance(atoms_num, base, pre_base, file_start, frames_num, each, start, end, \
             atom_type_1, atom_type_2, a_vec, b_vec, c_vec, file_name, work_dir):

  '''
  distance : calculate distance between atom type 1 and atom type 2 over different frame.

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
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    distance : 3-d float list, dim = frames_num*(number of atom_1)*(number of atom_2)
    atom_id_1 : 1-d int list
      atom_id_1 contains atom id of atom_type_1.
    atom_id_2 : 1-d int list
      atom_id_2 contains atom id of atom_type_2
  '''

  new_file_name = ''.join((work_dir, '/center_box.xyz'))
  center.center(atoms_num, base, pre_base, frames_num, a_vec, b_vec, c_vec, 'center_box', 0, file_name, work_dir)

  atom_id_1 = []
  atom_id_2 = []
  for i in range(atoms_num):
    line_i = linecache.getline(new_file_name, base+i+1)
    line_i_split = list_dic_op.str_split(line_i, ' ')
    if ( line_i_split[0] == atom_type_1 ):
      atom_id_1.append(i+1)
    if ( line_i_split[0] == atom_type_2 ):
      atom_id_2.append(i+1)

  frame_num_stat = int((end-start)/each+1)

  distance = []
  for i in range(frame_num_stat):
    distance_i = []
    for j in range(atoms_num):
      line_j = linecache.getline(new_file_name, (int((start-file_start)/each)+i)*(atoms_num+base)+j+base+1)
      line_j_split = list_dic_op.str_split(line_j, ' ')
      if ( line_j_split[0] == atom_type_1 ):
        coord_1 = []
        coord_2 = []
        for k in range(atoms_num):
          line_k = linecache.getline(new_file_name, (int((start-file_start)/each)+i)*(atoms_num+base)+k+base+1)
          line_k_split = list_dic_op.str_split(line_k, ' ')
          if ( line_k_split[0] == atom_type_2 ):
            coord_1.append([float(line_j_split[1]),float(line_j_split[2]),float(line_j_split[3].strip('\n'))])
            coord_2.append([float(line_k_split[1]),float(line_k_split[2]),float(line_k_split[3].strip('\n'))])
        dist = geometry_mod.geometry.calculate_distance(np.asfortranarray(coord_1, dtype='float32'), \
                                                        np.asfortranarray(coord_2, dtype='float32'), \
                                                        np.asfortranarray(a_vec, dtype='float32'), \
                                                        np.asfortranarray(b_vec, dtype='float32'), \
                                                        np.asfortranarray(c_vec, dtype='float32'))
        distance_i.append(list(dist))
    distance.append(distance_i)

  return distance, atom_id_1, atom_id_2

def rdf(distance, a_vec, b_vec, c_vec, r_increment, work_dir):

  '''
  rdf : get rdf between atom_1 and atom_2

  Args :
    distance : 3-d float list, dim = frames_num*(number of atom_1)*(number of atom_2)
    a_vec : 1d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec : 1d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec : 1d float list, dim = 3
      c_vec is the cell vector c.
      Example : [0.0, 0.0, 12.42]
    r_increment: float
      r_increment is the increment of r.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  vec = a_vec+b_vec+c_vec
  r_max = np.sqrt(np.dot(vec, vec))/2.0
  #vol = 4.0/3.0*np.pi*r_max**3
  vol = np.dot(a_vec, np.cross(b_vec, c_vec))
  data_num = int(r_max/r_increment)

  rdf_value, integral_value = \
  geometry_mod.geometry.rdf(distance, r_increment, vol, data_num)

  rdf_file = ''.join((work_dir, '/rdf_integral.csv'))
  with open(rdf_file ,'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['distance', 'rdf', 'int'])
    for i in range(data_num-1):
      writer.writerow([r_increment*(i+1), rdf_value[i], integral_value[i]])

def rdf_run(rdf_param, work_dir):

  '''
  rdf_run : the kernel function to run rdf function. It will call rdf.

  Args :
    rdf_param : dictionary
      rdf_param contains keywords used in rdf functions.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'traj_file' in rdf_param.keys() ):
    traj_file = rdf_param['traj_file']
    if ( os.path.exists(traj_file) ):
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      print ('Cannot find %s file' % (traj_file))
      exit()
  else:
    print ('No trajectory file found, please choose the traj_file')
    exit()

  if ( 'atom_type_pair' in rdf_param.keys() ):
    atom_type_pair = rdf_param['atom_type_pair']
    atom_1 = atom_type_pair[0]
    atom_2 = atom_type_pair[1]
  else:
    print ('No atom pair used to calculate bond length found, please set atom_pair')
    exit()

  if ( 'box' in rdf_param.keys() ):
    A_exist = 'A' in rdf_param['box'].keys()
    B_exist = 'B' in rdf_param['box'].keys()
    C_exist = 'C' in rdf_param['box'].keys()
  else:
    print ('No box information found, please set box')
    exit()

  if ( A_exist and B_exist and C_exist):
    box_A = rdf_param['box']['A']
    box_B = rdf_param['box']['B']
    box_C = rdf_param['box']['C']
  else:
    print ('Box information wrong, please check')
    exit()

  a_vec = [float(x) for x in box_A]
  b_vec = [float(x) for x in box_B]
  c_vec = [float(x) for x in box_C]

  if ( 'init_step' in rdf_param.keys() ):
    init_step = int(rdf_param['init_step'])
  else:
    init_step = start_frame_id

  if ( 'end_step' in rdf_param.keys() ):
    end_step = int(rdf_param['end_step'])
  else:
    end_step = start_frame_id + each

  if ( 'r_increment' in rdf_param.keys() ):
    r_increment = float(rdf_param['r_increment'])
  else:
    r_increment = 0.1

  dist, atom_1, atom_2 = distance(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, \
                                  end_step, atom_1, atom_2, a_vec, b_vec, c_vec, traj_file, work_dir)

  rdf(dist, a_vec, b_vec, c_vec, r_increment, work_dir)
