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
from CP2K_kit.analyze import check_analyze

def angle(atoms_num, pre_base_block, end_base_block, pre_base, start_frame_id, frames_num, each, \
          init_step, end_step, atom_type_1, atom_type_2, atom_type_3, traj_coord_file, work_dir):

  '''
  distance: calculate distance between atom type 1 and atom type 2 over different frame.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    pre_base_block: int
      pre_base_block is the number of lines before structure in a structure block.
    end_base_block: int
      end_base_block is the number of lines after structure in a structure block.
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

  atom_id_1 = []
  atom_id_2 = []
  atom_id_3 = []
  for i in range(atoms_num):
    line_i = linecache.getline(traj_coord_file, pre_base_block+pre_base+i+1)
    line_i_split = data_op.split_str(line_i, ' ')
    if ( line_i_split[0] == atom_type_1 ):
      atom_id_1.append(i+1)
    if ( line_i_split[0] == atom_type_2 ):
      atom_id_2.append(i+1)
    if ( line_i_split[0] == atom_type_3 ):
      atom_id_3.append(i+1)

  frame_num_stat = int((end_step-init_step)/each+1)

  angle = []
  for i in range(frame_num_stat):
    angle_i = []
    for j in atom_id_1:
      for k in atom_id_2:
        for l in atom_id_3:
          if ( l > j ):
            line_j_num = (int((init_step-start_frame_id)/each)+i)*(pre_base_block+atoms_num+end_base_block)+j+pre_base_block+pre_base
            line_j = linecache.getline(traj_coord_file, line_j_num)
            line_j_split = data_op.split_str(line_j, ' ', '\n')
            line_k_num = (int((init_step-start_frame_id)/each)+i)*(pre_base_block+atoms_num+end_base_block)+k+pre_base_block+pre_base
            line_k = linecache.getline(traj_coord_file, line_k_num)
            line_k_split = data_op.split_str(line_k, ' ', '\n')
            line_l_num = (int((init_step-start_frame_id)/each)+i)*(pre_base_block+atoms_num+end_base_block)+l+pre_base_block+pre_base
            line_l = linecache.getline(traj_coord_file, line_l_num)
            line_l_split = data_op.split_str(line_l, ' ', '\n')
            coord_1 = []
            coord_2 = []
            coord_3 = []
            coord_1.append([float(line_j_split[1]),float(line_j_split[2]),float(line_j_split[3])])
            coord_2.append([float(line_k_split[1]),float(line_k_split[2]),float(line_k_split[3])])
            coord_3.append([float(line_l_split[1]),float(line_l_split[2]),float(line_l_split[3])])
            ang = geometry_mod.geometry.calculate_angle(np.asfortranarray(coord_1, dtype='float32'), \
                                                        np.asfortranarray(coord_2, dtype='float32'), \
                                                        np.asfortranarray(coord_3, dtype='float32'))
            angle_i.append(ang[0])
    angle.append(angle_i)

  linecache.clearcache()

  return angle

def adf(angle, a_increment, work_dir):

  '''
  rdf: get rdf between atom type 1 and atom type 2

  Args:
    distance: 3-d float list, dim = frames_num*(number of atom_1)*(number of atom_2)
    a_increment: float
      a_increment is the increment of angle.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    adf_file: string
      adf_file contains adf information.
  '''

  a_max = 180.0
  data_num = int(a_max/a_increment)

  rho_a = geometry_mod.geometry.adf(angle, a_increment, data_num)

  adf_file = ''.join((work_dir, '/adf.csv'))
  with open(adf_file ,'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['angle(degree)', 'adf'])
    for i in range(data_num-1):
      writer.writerow([a_increment*(i+1), float(rho_a[i])/(len(angle)*len(angle[0]))])

  return adf_file

def adf_run(adf_param, work_dir):

  '''
  adf_run: the kernel function to run adf function. It will call adf.

  Args:
    adf_param: dictionary
      adf_param contains keywords used in adf functions.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  adf_param = check_analyze.check_adf_inp(adf_param)

  traj_coord_file = adf_param['traj_coord_file']
  atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
  traj_info.get_traj_info(traj_coord_file, 'coord_xyz')

  log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

  atom_type_pair = adf_param['atom_type_pair']
  atom_1 = atom_type_pair[0]
  atom_2 = atom_type_pair[1]
  atom_3 = atom_type_pair[2]

  init_step = adf_param['init_step']
  end_step = adf_param['end_step']
  a_increment = adf_param['a_increment']

  print ('ADF'.center(80, '*'), flush=True)
  print ('Analyze distribution of angular between %s and %s and %s' %(atom_1, atom_2, atom_3), flush=True)
  angle_whole = angle(atoms_num, pre_base_block, end_base_block, pre_base, start_frame_id, frames_num, each, \
                      init_step, end_step, atom_1, atom_2, atom_3,  traj_coord_file, work_dir)

  adf_file = adf(angle_whole, a_increment, work_dir)
  str_print = 'The adf file is written in %s' %(adf_file)
  print (data_op.str_wrap(str_print, 80), flush=True)
