#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import list_dic_op
from CP2K_kit.lib import rmsd_mod
from CP2K_kit.lib import statistic_mod

def rmsd(atoms_num, base, pre_base, each, atom_id, file_start, ref_frame, comp_frame_list, file_name):

  #Reference literature: J. Comput. Chem. 2004, 25, 1849-1857.

  '''
  rmsd : get rmsd of choosed atoms in md trajectory.

  Args :
    atoms_num : int
      atoms_num is the number of atoms of the system.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    each : int
      each is printing frequency of md.
    atom_id : int list
      atom_id is the id of atoms to be analyzed.
      Example : [1,2,3,7,8]
    file_start : int
      file_start is the starting frame in trajectory file.
    ref_frame : int
      ref_frame is the reference frame.
    comp_frame_list : 1-d int list
      comp_frame_list is comparring frames.
    file_name : string
      file_name is the name of trajectory file used to analyze.
  Returns :
    rmsd_value_list : 1-d float list
  '''

  coord_ref = np.asfortranarray(np.zeros((len(atom_id),3)),dtype='float32')
  coord_comp = np.asfortranarray(np.zeros((len(atom_id),3)),dtype='float32')

  for i in range(len(atom_id)):
    line_i = linecache.getline(file_name, int((ref_frame-file_start)/each)*(atoms_num+base)+atom_id[i]+base+pre_base)
    line_i_split = list_dic_op.str_split(line_i, ' ')
    coord_ref[i,0] = float(line_i_split[1])
    coord_ref[i,1] = float(line_i_split[2])
    coord_ref[i,2] = float(line_i_split[3].strip('\n'))

  coord_ref_center = np.asfortranarray(np.zeros(3),dtype='float32')
  for i in range(3):
    coord_ref_center[i] = statistic_mod.statistic.numerical_average(coord_ref[:,i],len(atom_id))

  rmsd_value_list = []

  for m in range(len(comp_frame_list)):
    for i in range(len(atom_id)):
      line_mi = linecache.getline(file_name, int((comp_frame_list[m]-file_start)/each)*(atoms_num+base)+atom_id[i]+base+pre_base)
      line_mi_split = list_dic_op.str_split(line_mi, ' ')
      coord_comp[i,0] = float(line_mi_split[1])
      coord_comp[i,1] = float(line_mi_split[2])
      coord_comp[i,2] = float(line_mi_split[3].strip('\n'))

    coord_comp_center = np.asfortranarray(np.zeros(3),dtype='float32')

    for i in range(3):
      coord_comp_center[i] = statistic_mod.statistic.numerical_average(coord_comp[:,i],len(atom_id))

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

  return rmsd_value_list

def rmsd_run(rmsd_param, work_dir):

  '''
  rmsd_run : the kernel function to run rmsd function.

  Args :
    rmsd_param : dictionary
      rmsd_param contains keywords used in rmsd functions.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'traj_file' in rmsd_param.keys() ):
    traj_file = rmsd_param['traj_file']
    if ( os.path.exists(traj_file) ):
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      log_info.log_error('%s file does not exist' %(traj_file))
  else:
    log_info.log_error('No trajectory found, please set analyze/rmsd/traj_file')
    exit()

  if ( 'atom_id' in rmsd_param.keys() ):
    atom_id = rmsd_param['atom_id']
    atom_id_list = list_dic_op.get_id_list(atom_id)
  else:
    log_info.log_error('No atom id found, please set analyze/rmsd/atom_id')
    exit()

  if ( 'ref_frame' in rmsd_param.keys() ):
    ref_frame = int(rmsd_param['ref_frame'])
  else:
    log_info.log_error('No reference frame found, please set analyze/rmsd/ref_frame')
    exit()

  if ( 'compare_frame' in rmsd_param.keys() ):
    compare_frame = rmsd_param['compare_frame']
    compare_frame_list = list_dic_op.get_id_list(compare_frame)
  else:
    log_info.log_error('No compare frame found, please set analyze/rmsd/compare_frame')
    exit()

  rmsd_value = rmsd(atoms_num, base, pre_base, each, atom_id_list, \
                    start_frame_id, ref_frame, compare_frame_list, traj_file)

  rmsd_file = ''.join((work_dir, '/rmsd.csv'))
  with open(rmsd_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'rmsd'])
    for i in range(len(rmsd_value)):
      writer.writerow([i*time_step*each, rmsd_value[i]])
