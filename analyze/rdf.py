#!/usr/bin/env python

import os
import csv
import math
import linecache
import numpy as np
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import call
from CP2K_kit.lib import geometry_mod
from CP2K_kit.analyze import center
from CP2K_kit.analyze import check_analyze

def distance(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, end_step, \
             atom_type_1, atom_type_2, a_vec, b_vec, c_vec, traj_coord_file, work_dir):

  '''
  distance: calculate distance between atom type 1 and atom type 2 over different frame.

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
      Example : [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    traj_coord_file: string
      file_name is the name of coordination trajectory file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    distance: 3-d float list, dim = frames_num*(number of atom_1)*(number of atom_2)
              if atom_type_1 = atom_type_2, dim = frames_num*(number of atom_1 - 1)*(number of atom_2 - 1)
    atom_id_1: 1-d int list
      atom_id_1 contains atom id of atom_type_1.
    atom_id_2: 1-d int list
      atom_id_2 contains atom id of atom_type_2
  '''

  center_file = center.center(atoms_num, base, pre_base, frames_num, a_vec,  b_vec, \
                c_vec, 'center_box', 0, traj_coord_file, work_dir, 'center.xyz')

  atom_id_1 = []
  atom_id_2 = []
  for i in range(atoms_num):
    line_i = linecache.getline(center_file, base+i+1)
    line_i_split = data_op.str_split(line_i, ' ')
    if ( line_i_split[0] == atom_type_1 ):
      atom_id_1.append(i+1)
    if ( line_i_split[0] == atom_type_2 ):
      atom_id_2.append(i+1)

  frame_num_stat = int((end_step-init_step)/each+1)

  distance = []
  for i in range(frame_num_stat):
    distance_i = []
    for j in range(atoms_num):
      line_j = linecache.getline(center_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+base)+j+base+1)
      line_j_split = data_op.str_split(line_j, ' ')
      if ( line_j_split[0] == atom_type_1 ):
        coord_1 = []
        coord_2 = []
        for k in range(atoms_num):
          line_k = linecache.getline(center_file, (int((init_step-start_frame_id)/each)+i)*(atoms_num+base)+k+base+1)
          line_k_split = data_op.str_split(line_k, ' ')
          if ( line_k_split[0] == atom_type_2 and j != k ):
            coord_1.append([float(line_j_split[1]),float(line_j_split[2]),float(line_j_split[3].strip('\n'))])
            coord_2.append([float(line_k_split[1]),float(line_k_split[2]),float(line_k_split[3].strip('\n'))])
        dist = geometry_mod.geometry.calculate_distance(np.asfortranarray(coord_1, dtype='float32'), \
                                                        np.asfortranarray(coord_2, dtype='float32'), \
                                                        np.asfortranarray(a_vec, dtype='float32'), \
                                                        np.asfortranarray(b_vec, dtype='float32'), \
                                                        np.asfortranarray(c_vec, dtype='float32'))
        distance_i.append(list(dist))
    distance.append(distance_i)

  cmd = 'rm %s' %(center_file)
  call.call_simple_shell(work_dir, cmd)

  return distance, atom_id_1, atom_id_2

def rdf(distance, a_vec, b_vec, c_vec, r_increment, work_dir):

  '''
  rdf: get rdf between atom type 1 and atom type 2

  Args:
    distance: 3-d float list, dim = frames_num*(number of atom_1)*(number of atom_2)
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example: [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Exampl : [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    r_increment: float
      r_increment is the increment of r.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    rdf_file: string
      rdf_file contains rdf information.
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
    writer.writerow(['distance(Ang)', 'rdf', 'int'])
    for i in range(data_num-1):
      writer.writerow([r_increment*(i+1), rdf_value[i], integral_value[i]])

  return rdf_file

def rdf_run(rdf_param, work_dir):

  '''
  rdf_run: the kernel function to run rdf function. It will call rdf.

  Args:
    rdf_param: dictionary
      rdf_param contains keywords used in rdf functions.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  rdf_param = check_analyze.check_rdf_inp(rdf_param)

  traj_coord_file = rdf_param['traj_coord_file']
  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
  traj_info.get_traj_info(traj_coord_file, 'coord')

  log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

  atom_type_pair = rdf_param['atom_type_pair']
  atom_1 = atom_type_pair[0]
  atom_2 = atom_type_pair[1]

  a_vec = rdf_param['box']['A']
  b_vec = rdf_param['box']['B']
  c_vec = rdf_param['box']['C']

  init_step = rdf_param['init_step']
  end_step = rdf_param['end_step']
  r_increment = rdf_param['r_increment']

  print ('RDF'.center(80, '*'), flush=True)
  print ('Analyze radial distribution function between %s and %s' %(atom_1, atom_2), flush=True)
  dist, atom_id_1, atom_id_2 = distance(atoms_num, base, pre_base, start_frame_id, frames_num, each, init_step, \
                                        end_step, atom_1, atom_2, a_vec, b_vec, c_vec, traj_coord_file, work_dir)

  rdf_file = rdf(dist, a_vec, b_vec, c_vec, r_increment, work_dir)
  str_print = 'The rdf file is written in %s' %(rdf_file)
  print (data_op.str_wrap(str_print, 80), flush=True)
