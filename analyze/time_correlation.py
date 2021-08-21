#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from CP2K_kit.tools import atom
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.lib import statistic_mod
from CP2K_kit.lib import dynamic_mod

#The unit of velocity in CP2K trajectory is Bohr/au_t
#We transfer velocity in cm/s
#The unit of intensity is cm^2/s

def time_corr_func(atoms_num, base, pre_base, each, file_start, time_step, start, end, \
                   time_num, atom_id, file_name, work_dir, normalize=1, return_func=False):

  '''
  time_corr_func : calculate time correlation function.

  Args :
    atoms_num : int
      atoms_num is the number of atoms in trajectory file.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    each : int
      each is printing frequency of md.
    file_start : int
      file_start is the starting frame in trajectory file.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    time_num : int
      time_num is the max correlation frame number.
    atom_id : int list
      atom_id is the id of atoms to be analyzed.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    normalize : int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    return_func : bool
      return_func is whether to return time correlation function.
  Returns :
    data_tcf : 1-d float array
      data_tcf is the time correlation function for a physical quantity.
  '''

  frame_num_stat = int((end-start)/each+1)
  data = np.asfortranarray(np.zeros((frame_num_stat,len(atom_id),3)),dtype='float32')
  for i in range(frame_num_stat):
    for j in range(len(atom_id)):
      line_ij = linecache.getline(file_name,(atoms_num+base)*(int((start-file_start)/each)+i)+base+atom_id[j])
      line_ij_split = data_op.str_split(line_ij, ' ')
      data[i,j,0] = float(line_ij_split[1])
      data[i,j,1] = float(line_ij_split[2])
      data[i,j,2] = float(line_ij_split[3].strip('\n'))
  data_tcf = dynamic_mod.dynamic.time_correlation(data, time_num, normalize)

  tcf_file = ''.join((work_dir, '/tcf.csv'))
  with open(tcf_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'acf'])
    for i in range(len(data_tcf)):
      writer.writerow([i*time_step,data_tcf[i]])

  if return_func:
    return data_tcf

def time_corr_mode_func(atoms_num, base, pre_base, each, file_start, time_step, start, end, time_num, cluster_group_id,\
                         pos_file, vel_file, a_vec, b_vec, c_vec, work_dir, normalize=1, return_func=False):

  '''
  time_corr_mode_func : calculate mode time correlation function.

  Args :
    atoms_num : int
      atoms_num is the number of atoms in trajectory file.
    base : int
      base is the number of lines before structure in a structure block.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    each : int
      each is printing frequency of md.
    file_start : int
      file_start is the starting frame in trajectory file.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    time_num : int
      time_num is the max correlation frame number.
    cluster_group_id : 3-d int list
      cluster_group_id is the id of atoms in the molecule group.
    pos_file : string
      pos_file is the position trajectory file.
    vel_file : string
      vel_file is the velocity trajectory file.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    normalize : int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    return_func : bool
      return_func is whether to return time correlation function.
  Returns :
    data_q1_tcf : 1-d float array
      data_q1_tcf is the time correlation function for a physical quantity of q1 mode.
    data_q2_tcf : 1-d float array
      data_q2_tcf is the time correlation function for a physical quantity of q2 mode.
    data_q3_tcf : 1-d float array
      data_q3_tcf is the time correlation function for a physical quantity of q3 mode.
  '''

  frame_num_stat = int((end-start)/each+1)

  Q1_data = np.asfortranarray(np.zeros((frame_num_stat,len(cluster_group_id[0]),3)),dtype='float32')
  Q2_data = np.asfortranarray(np.zeros((frame_num_stat,len(cluster_group_id[0]),3)),dtype='float32')
  Q3_data = np.asfortranarray(np.zeros((frame_num_stat,len(cluster_group_id[0]),3)),dtype='float32')

  for i in range(frame_num_stat):
    for j in range(len(cluster_group_id[i])):
      pos_data = np.asfortranarray(np.zeros((len(cluster_group_id[i][j]),3)),dtype='float32')
      vel_data = np.asfortranarray(np.zeros((len(cluster_group_id[i][j]),3)),dtype='float32')
      element = []
      for k in range(len(cluster_group_id[i][j])):
        #Dump coordinate
        line_k = linecache.getline(pos_file, (atoms_num+base)*(int((start-file_start)/each)+i)+cluster_group_id[i][j][k]+base+pre_base)
        line_k_split = data_op.str_split(line_k, ' ')
        element.append(line_k_split[0])
        pos_data[k,0] = float(line_k_split[1])
        pos_data[k,1] = float(line_k_split[2])
        pos_data[k,2] = float(line_k_split[3].strip('\n'))

        #Dump velocity
        line_k = linecache.getline(vel_file, (atoms_num+base)*(int((start-file_start)/each)+i)+cluster_group_id[i][j][k]+base+pre_base)
        line_k_split = data_op.str_split(line_k, ' ')
        vel_data[k,0] = float(line_k_split[1])
        vel_data[k,1] = float(line_k_split[2])
        vel_data[k,2] = float(line_k_split[3].strip('\n'))

      atom_number, atom_mass = atom.get_atom_mass(element)
      mass_array = np.asfortranarray(atom_mass, dtype='float32')
      q1, q2, q3 = statistic_mod.statistic.data_mode(pos_data,vel_data,a_vec,b_vec,c_vec,mass_array)
      Q1_data[i,j,0] = q1[0]
      Q1_data[i,j,1] = q1[1]
      Q1_data[i,j,2] = q1[2]
      Q2_data[i,j,0] = q2[0]
      Q2_data[i,j,1] = q2[1]
      Q2_data[i,j,2] = q2[2]
      Q3_data[i,j,0] = q3[0]
      Q3_data[i,j,1] = q3[1]
      Q3_data[i,j,2] = q3[2]

  data_q1_tcf = dynamic_mod.dynamic.time_correlation(Q1_data, time_num, normalize)
  data_q2_tcf = dynamic_mod.dynamic.time_correlation(Q2_data, time_num, normalize)
  data_q3_tcf = dynamic_mod.dynamic.time_correlation(Q3_data, time_num, normalize)

  tcf_q1_file = ''.join((work_dir, '/tcf_q1.csv'))
  with open(tcf_q1_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'acf'])
    for i in range(len(data_q1_tcf)):
      writer.writerow([i*time_step,data_q1_tcf[i]])

  tcf_q2_file = ''.join((work_dir, '/tcf_q2.csv'))
  with open(tcf_q2_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'acf'])
    for i in range(len(data_q2_tcf)):
      writer.writerow([i*time_step,data_q2_tcf[i]])

  tcf_q3_file = ''.join((work_dir, '/tcf_q3.csv'))
  with open(tcf_q3_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'acf'])
    for i in range(len(data_q3_tcf)):
      writer.writerow([i*time_step,data_q3_tcf[i]])

  if return_func:
    return data_q1_tcf, data_q2_tcf, data_q3_tcf

def time_corr_run(time_corr_param, work_dir):

  '''
  time_corr_run : the kernel function to run time correlation function.

  Args :
    time_corr_param : dictionary
      time_corr_param contains keywords used in time correlation functions.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'traj_file' in time_corr_param.keys() ):
    traj_file = time_corr_param['traj_file']
    if ( os.path.exists(traj_file) ):
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      log_info.log_error('%s file does not exists' %(traj_file))
      exit()
  else:
    log_info.log_error('No trajectory file found, please set analyze/time_correlation/traj_file')
    exit()

  if ( 'atom_id' in time_corr_param.keys() ):
    atom_id = time_corr_param['atom_id']
    atom_id_list = data_op.get_id_list(atom_id)
  else:
    log_info.log_error('No atom id found, please set analyze/time_correlation/atom_id')
    exit()

  if ( 'max_frame_corr' in time_corr_param.keys() ):
    max_frame_corr = int(time_corr_param['max_frame_corr'])
  else:
    max_frame_corr = int(frames_num/3)

  if ( 'init_step' in time_corr_param.keys() ):
    init_step = int(time_corr_param['init_step'])
  else:
    init_step = start_frame_id

  if ( 'end_step' in time_corr_param.keys() ):
    end_step = int(time_corr_param['end_step'])
  else:
    end_step = start_frame_id + each

  if ( 'normalize' in time_corr_param.keys() ):
    normalize = int(time_corr_param['normalize'])
  else:
    normalize = 1

  time_corr_func(atoms_num, base, pre_base, each, start_frame_id, time_step, init_step, \
                 end_step, max_frame_corr, atom_id_list, traj_file, work_dir, normalize)
