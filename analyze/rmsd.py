#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import data_op
from CP2K_kit.analyze import check_analyze
from CP2K_kit.lib import rmsd_mod
from CP2K_kit.lib import statistic_mod

def rmsd(atoms_num, base, pre_base, each, atom_id, start_frame_id, ref_frame, comp_frame_list, traj_coord_file):

  #Reference literature: J. Comput. Chem. 2004, 25, 1849-1857.

  '''
  rmsd: get rmsd of choosed atoms in md trajectory.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    each: int
      each is printing frequency of md.
    atom_id: int list
      atom_id is the id of atoms.
      Example: [1,2,3,7,8]
    start_frame_id: int
      start_frame_id is the starting frame id in trajectory file.
    ref_frame: int
      ref_frame is the reference frame.
    comp_frame_list: 1-d int list
      comp_frame_list is comparring frames.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
  Returns:
    rmsd_value_list: 1-d float list
      rmsd_value_list is the list of rmsd value.
  '''

  coord_ref = np.asfortranarray(np.zeros((len(atom_id),3)),dtype='float32')
  coord_comp = np.asfortranarray(np.zeros((len(atom_id),3)),dtype='float32')

  for i in range(len(atom_id)):
    line_i = linecache.getline(traj_coord_file, int((ref_frame-start_frame_id)/each)*(atoms_num+base)+atom_id[i]+base+pre_base)
    line_i_split = data_op.str_split(line_i, ' ')
    coord_ref[i,0] = float(line_i_split[1])
    coord_ref[i,1] = float(line_i_split[2])
    coord_ref[i,2] = float(line_i_split[3].strip('\n'))

  coord_ref_center = np.asfortranarray(np.zeros(3),dtype='float32')
  for i in range(3):
    value_avg, sigma = statistic_mod.statistic.numerical_average(coord_ref[:,i],len(atom_id))
    coord_ref_center[i] = value_avg

  rmsd_value_list = []

  for m in range(len(comp_frame_list)):
    for i in range(len(atom_id)):
      line_mi = linecache.getline(traj_coord_file, int((comp_frame_list[m]-start_frame_id)/each)*(atoms_num+base)+atom_id[i]+base+pre_base)
      line_mi_split = data_op.str_split(line_mi, ' ')
      coord_comp[i,0] = float(line_mi_split[1])
      coord_comp[i,1] = float(line_mi_split[2])
      coord_comp[i,2] = float(line_mi_split[3].strip('\n'))

    coord_comp_center = np.asfortranarray(np.zeros(3),dtype='float32')

    for i in range(3):
      value_avg, sigma = statistic_mod.statistic.numerical_average(coord_comp[:,i],len(atom_id))
      coord_comp_center[i] = value_avg

    cov_matrix = rmsd_mod.rmsd.get_cov_matrix(coord_comp,coord_ref,coord_comp_center,coord_ref_center)
    quart_matrix = rmsd_mod.rmsd.quarternion_rotate(cov_matrix)

    eigen = np.linalg.eig(quart_matrix)
    eigen_value = eigen[0]
    eigen_vector = eigen[1]
    eigen_max = max(eigen_value)
    max_index = list(eigen_value).index(max(eigen_value))
    quart_vec = eigen_vector[max_index]

    #If you need rotate matrix, you could print it or return it

    rotate_max = rmsd_mod.rmsd.quart_to_rot(quart_vec)
    rmsd_value = rmsd_mod.rmsd.get_rmsd(coord_comp,coord_ref,coord_comp_center,coord_ref_center,eigen_max)
    rmsd_value_list.append(rmsd_value)

  linecache.clearcache()

  return rmsd_value_list

def rmsd_run(rmsd_param, work_dir):

  '''
  rmsd_run: the kernel function to run rmsd function.

  Args:
    rmsd_param: dictionary
      rmsd_param contains keywords used in rmsd functions.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  rmsd_param = check_analyze.check_rmsd_inp(rmsd_param)

  traj_coord_file = rmsd_param['traj_coord_file']
  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
  traj_info.get_traj_info(traj_coord_file, 'coord')

  log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

  atom_id = rmsd_param['atom_id']
  ref_frame = rmsd_param['ref_frame']
  compare_frame = rmsd_param['compare_frame']

  print ('RMSD'.center(80, '*'), flush=True)
  print ('Calculate root mean square deviation based on reference frame %d' %(ref_frame), flush=True)

  rmsd_value = rmsd(atoms_num, base, pre_base, each, atom_id, start_frame_id, ref_frame, compare_frame, traj_coord_file)

  rmsd_file = ''.join((work_dir, '/rmsd.csv'))
  with open(rmsd_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'rmsd'])
    for i in range(len(rmsd_value)):
      writer.writerow([i*time_step*each, rmsd_value[i]])

  str_print = 'The rmsd vs time is written in %s' %(rmsd_file)
  print (data_op.str_wrap(str_print, 80), flush=True)

