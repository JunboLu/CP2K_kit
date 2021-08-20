#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from CP2K_kit.tools import log_info
from CP2K_kit.tools import list_dic_op
from CP2K_kit.tools import numeric
from CP2K_kit.tools import traj_info
from CP2K_kit.analyze import time_correlation
from CP2K_kit.analyze import geometry
from CP2K_kit.lib import statistic_mod

#The unit of velocity is CP2K trajectory is Bohr/au_t

def power_spectrum(atoms_num, base, pre_base, each, file_start, time_step, start, end, time_num,
                   lower_wave, upper_wave, increment_wave, atom_id, file_name, normalize, work_dir):

  '''
  power_spectrum : calculate power_spectrum of the system

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
    lower_wave : float
      lower_wave is the starting wave value. The unit is cm^-1.
    upper_wave : float
      upper_wave is the ending wave value. The unit is cm^-1.
    increment_wave : float
      increment_wave is the increment of wave. The unit is cm^-1.
    atom_id : int list
      atom_id is the id of atoms to be analyzed.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    normalize : int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  vel_tcf = time_correlation.time_corr_func(atoms_num, base, pre_base, each, file_start, time_step, start, end,
                                            time_num, atom_id, file_name, work_dir, normalize, return_func=True)

  data_num = int((upper_wave-lower_wave)/increment_wave+1)
  wave_num = np.asfortranarray(np.zeros(data_num),dtype='float32')
  time = np.asfortranarray(np.zeros(len(vel_tcf)),dtype='float32')

  for i in range(data_num):
    wave_num[i] = lower_wave+increment_wave*i

  for i in range(len(time)):
    time[i]=i*time_step*each

  if (normalize == 0):
    time_step = time_step*1E-15

  intensity = \
  statistic_mod.statistic.fourier_transform(vel_tcf, time, wave_num, time_step)

  intensity_fit = numeric.savitzky_golay(intensity,201,3)

  freq_int_file = ''.join((work_dir, '/freq_intensity.csv'))
  with open(freq_int_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['wave', 'intensity', 'intensity_fit'])
    for i in range(data_num):
      writer.writerow([wave_num[i], intensity[i], intensity_fit[i]])

def power_spectrum_mode(atoms_num, base, pre_base, each, file_start, time_step, start, end, time_num, cluster_group_id,
                        lower_wave, upper_wave, increment_wave, pos_file, vel_file, a_vec, b_vec, c_vec, normalize, work_dir):

  #Reference literature: Chem. Phys. 1986, 106, 205-212.

  '''
  power_spectrum : calculate mode power_spectrum of the system

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
    cluster_group_id : 1-d int list
      cluster_group_id is the id of first atoms in the molecules in the group.
    lower_wave : float
      lower_wave is the starting wave value. The unit is cm^-1.
    upper_wave : float
      upper_wave is the ending wave value. The unit is cm^-1.
    increment_wave : float
      increment_wave is the increment of wave. The unit is cm^-1.
    pos_file : string
      pos_file is the position trajectory file.
    vel_file : string
      vel_file is the velocity trajectory file.
    normalize : int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  Q1_vel_tcf, Q2_vel_tcf, Q3_vel_tcf = \
  time_correlation.time_corr_mode_func(atoms_num, base, pre_base, each, file_start, time_step, start, end, time_num, cluster_group_id, \
                                       pos_file, vel_file, a_vec, b_vec, c_vec, work_dir, normalize, return_func=True)

  data_num = int((upper_wave-lower_wave)/increment_wave+1)
  wave_num = np.asfortranarray(np.zeros(data_num),dtype='float32')
  time = np.asfortranarray(np.zeros(len(Q1_vel_tcf)),dtype='float32')

  for i in range(data_num):
    wave_num[i] = lower_wave+increment_wave*i

  for i in range(len(time)):
    time[i]=i*time_step

  if (normalize == 0):
    time_step = time_step*1E-15

  Q1_intensity = \
  statistic_mod.statistic.fourier_transform(Q1_vel_tcf, time, wave_num, time_step)
  Q1_intensity_fit = numeric.savitzky_golay(Q1_intensity,201,3)

  freq_int_q1_file = ''.join((work_dir, '/freq_intensity_q1.csv'))
  with open(freq_int_q1_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['wave', 'intensity', 'intensity_fit'])
    for i in range(data_num):
      writer.writerow([wave_num[i], Q1_intensity[i], Q1_intensity_fit[i]])

  Q2_intensity = \
  statistic_mod.statistic.fourier_transform(Q2_vel_tcf, time, wave_num, time_step)
  Q2_intensity_fit = numeric.savitzky_golay(Q2_intensity,201,3)

  freq_int_q2_file = ''.join((work_dir, '/freq_intensity_q2.csv'))
  with open(freq_int_q2_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['wave', 'intensity', 'intensity_fit'])
    for i in range(data_num):
      writer.writerow([wave_num[i], Q2_intensity[i], Q2_intensity_fit[i]])

  Q3_intensity = \
  statistic_mod.statistic.fourier_transform(Q3_vel_tcf, time, wave_num, time_step)
  Q3_intensity_fit = numeric.savitzky_golay(Q3_intensity,201,3)

  freq_int_q3_file = ''.join((work_dir, '/freq_intensity_q3.csv'))
  with open(freq_int_q3_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['wave', 'intensity', 'intensity_fit'])
    for i in range(data_num):
      writer.writerow([wave_num[i], Q3_intensity[i], Q3_intensity_fit[i]])

def power_spectrum_run(spectrum_param, work_dir):

  '''
  power_spectrum_run : the kernel function to run power_spectrum function.

  Args :
    spectrum_param : dictionary
      spectrum_param contains keywords used in power_spectrum functions.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'type' in spectrum_param.keys() ):
    spec_type = spectrum_param['type']
  else:
    log_info('No type found, please set analyze/power_spectrum/type')
    exit()

  if ( spec_type == 'water_mode' or spec_type == 'hydration_mode' ):
    if ( 'traj_pos_file' in spectrum_param.keys() ):
      traj_pos_file = spectrum_param['traj_pos_file']
      if ( os.path.exists(traj_pos_file) ):
        pass
      else:
        log_info.log_error('%s file does not exist' %(traj_pos_file))
        exit()
    else:
      log_info.log_error('No position trajectory found, please set analyze/power_spectrum/traj_pos_file')
      exit()

  if ( 'traj_vel_file' in spectrum_param.keys() ):
    traj_vel_file = spectrum_param['traj_vel_file']
    if ( os.path.exists(traj_vel_file) ):
      pass
    else:
      log_info.log_error('%s file does not exist' %(traj_vel_file))
  else:
    log_info.log_error('No velocity trajectory found, please set analyze/power_spectrum/traj_vel_file')

  if ( 'init_step' in spectrum_param.keys() ):
    init_step = int(spectrum_param['init_step'])
  else:
    log_info.log_error('No inititial step found, please set analyze/power_spectrum/init_step')
    exit()

  if ( 'end_step' in spectrum_param.keys() ):
    end_step = int(spectrum_param['end_step'])
  else:
    log_info.log_error('No end step found, please set analyze/power_spectrum/end_step')
    exit()

  if ( 'max_frame_corr' in spectrum_param.keys() ):
    max_frame_corr = int(spectrum_param['max_frame_corr'])
  else:
    max_frame_corr = int(frames_num/3)

  if ( 'start_wave' in spectrum_param.keys() ):
    start_wave = float(spectrum_param['start_wave'])
  else:
    start_wave = 0.0

  if ( 'end_wave' in spectrum_param.keys() ):
    end_wave = float(spectrum_param['end_wave'])
  else:
    end_wave = 4000.0

  if ( 'increment_wave' in spectrum_param.keys() ):
    increment_wave = float(spectrum_param['increment_wave'])
  else:
    increment_wave = 1.0

  if ( 'normalize' in spectrum_param.keys() ):
    normalize = int(spectrum_param['normalize'])
  else:
    normalize = 1

  if ( spec_type == 'general' ):
    if ( 'atom_id' in spectrum_param.keys() ):
      atom_id_list = list_dic_op.get_id_list(spectrum_param['atom_id'])
    else:
      log_info.log_error('No atom id found, please set analyze/power_spectrum/atom_id')
      exit()

    atoms_num, base, pre_base, frame_num, each, start_id, end_id, time_step = \
    traj_info.get_traj_info(traj_vel_file)

    power_spectrum(atoms_num, base, pre_base, each, start_id, time_step, init_step, end_step, max_frame_corr, \
                   start_wave, end_wave, increment_wave, atom_id_list, traj_vel_file, normalize, work_dir)

  else:
    atoms_num_p, base_p, pre_base_p, frame_num_p, each_p, start_id_p, end_id_p, time_step_p = \
    traj_info.get_traj_info(traj_pos_file)
    atoms_num_v, base_v, pre_base_v, frame_num_v, each_v, start_id_v, end_id_v, time_step_v = \
    traj_info.get_traj_info(traj_vel_file)

    if ( 'box' in spectrum_param.keys() ):
      A_exist = 'A' in spectrum_param['box'].keys()
      B_exist = 'B' in spectrum_param['box'].keys()
      C_exist = 'C' in spectrum_param['box'].keys()
    else:
      log_info.log_error('No box found, please set analyze/power_spectrum/box')
      exit()

    if ( A_exist and B_exist and C_exist):
      box_A = spectrum_param['box']['A']
      box_B = spectrum_param['box']['B']
      box_C = spectrum_param['box']['C']
    else:
      log_info.log_error('Box setting error, please check analzye/power_spectrum/box')
      exit()

    a_vec = [float(x) for x in box_A]
    b_vec = [float(x) for x in box_B]
    c_vec = [float(x) for x in box_C]

    if ( spec_type == 'water_mode' ):
      if (start_id_v > start_id_p):
        atoms_num, base, pre_base, frame_num, each, start_id, end_id, \
        time_step, exclude_group_id, group_atom_1_id, group_atoms_mass = \
        traj_info.get_traj_info(traj_vel_file, [['O','H','H']], True)
      else:
        atoms_num, base, pre_base, frame_num, each, start_id, end_id, \
        time_step, exclude_group_id, group_atom_1_id, group_atoms_mass = \
        traj_info.get_traj_info(traj_pos_file, [['O','H','H']], True)

      cluster_id = []
      for i in range(int(end_step)-int(init_step)+1):
        cluster_frame_id = []
        for j in range(len(group_atom_1_id[0])):
          cluster_frame_id.append([group_atom_1_id[0][j], group_atom_1_id[0][j]+1, group_atom_1_id[0][j]+2])
        cluster_id.append(cluster_frame_id)
      power_spectrum_mode(atoms_num, base, pre_base, each, start_id, time_step, init_step, end_step, max_frame_corr, cluster_id, \
                          start_wave, end_wave, increment_wave, traj_pos_file, traj_vel_file, a_vec, b_vec, c_vec, normalize, work_dir)

    elif ( spec_type == 'hydration_mode' ):
      if (start_id_v > start_id_p):
        atoms_num, base, pre_base, frame_num, each, start_id, end_id, time_step = \
        traj_info.get_traj_info(traj_vel_file)
      else:
        atoms_num, base, pre_base, frame_num, each, start_id, end_id, time_step = \
        traj_info.get_traj_info(traj_pos_file)

      if ( 'hyd_shell_dist' in spectrum_param.keys() ):
        hyd_shell_dist = float(spectrum_param['hyd_shell_dist'])
      else:
        log_info.log_error('No hydration shell distance found, please set analyze/power_spectrum/hyd_shell_num')
        exit()

      if ( 'dist_conv' in spectrum_param.keys() ):
        dist_conv = float(spectrum_param['dist_conv'])
      else:
        dist_conv = 0.3

      if ( 'atom_pair' in spectrum_param.keys() ):
        atom_pair = spectrum_param['atom_pair']
        atom_1 = atom_pair[0]
        atom_2 = atom_pair[1]
      else:
        log_info.log_error('No atom pair found, please set analyze/power_spectrum/atom_pair')
        exit()

      first_shell_id, dist = geometry.first_shell(atoms_num, base, pre_base, start_id, frames_num, each, init_step, end_step, \
                                                  atom_1, atom_2, a_vec, b_vec, c_vec, traj_pos_file, hyd_shell_dist, dist_conv)

      cluster_id = []
      for i in range(int(end_step)-int(init_step)+1):
        angle_order = geometry.order_angle(first_shell_id[i][0][0], first_shell_id[i][0][1:len(first_shell_id[i][0])-1], \
                                           i, init_step, traj_pos_file)
        angle_order_cyclic = angle_order+angle_order
        cluster_frame_id = []
        for j in range(len(angle_order)):
          group_id = []
          group_id.append(first_shell_id[i][0][0])
          group_id.append(angle_order_cyclic[j])
          group_id.append(angle_order_cyclic[j+1])
          cluster_frame_id.append(group_id)
        cluster_id.append(cluster_frame_id)

      power_spectrum_mode(atoms_num, base, pre_base, each, start_id, time_step, init_step, end_step, max_frame_corr, cluster_id, \
                          start_wave, end_wave, increment_wave, traj_pos_file, traj_vel_file, a_vec, b_vec, c_vec, normalize, work_dir)
