#!/usr/bin/env python

import os
import csv
import numpy as np
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import numeric
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import traj_tools
from CP2K_kit.analyze import check_analyze
from CP2K_kit.analyze import time_correlation
from CP2K_kit.analyze import geometry
from CP2K_kit.lib import statistic_mod

#The unit of velocity is CP2K trajectory is Bohr/au_t

def power_spectrum(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, end_step, \
                   max_frame_corr, lower_wave, upper_wave, increment_wave, atom_id, traj_vel_file, normalize, work_dir):

  '''
  power_spectrum: calculate power_spectrum of the system

  Args :
    atoms_num: int
      atoms_num is the number of atoms in trajectory file.
    pre_base_block: int
      pre_base_block is the number of lines before structure in a structure block.
    end_base_block: int
      end_base_block is the number of lines after structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    each: int
      each is printing frequency of md.
    start_frame_id: int
      start_frame_id is the starting frame id in the trajectory file.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    max_frame_corr: int
      max_frame_corr is the max number of correlation frames.
    lower_wave: float
      lower_wave is the starting wave value. The unit is cm^-1.
    upper_wave: float
      upper_wave is the ending wave value. The unit is cm^-1.
    increment_wave: float
      increment_wave is the increment of wave. The unit is cm^-1.
    atom_id: int list
      atom_id is the id of atoms.
    traj_vel_file: string
      traj_vel_file is the name of velocity trajectory file.
    normalize: int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns :
    wave_num: 1-d float list
      wave_num is the list of wave number.
    intensity: 1-d float array
      intensity is the list of intensity for each wave number.
    intensity_fit: 1-d float array
      intensity_fit is the list of fitted intensity for each wave number.
  '''

  print ('Calculate velocity-velocity auto correlation function at first', flush=True)
  vel_tcf, tcf_file = time_correlation.time_corr_func(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, \
                      time_step, init_step, end_step, max_frame_corr, atom_id, traj_vel_file, work_dir, 'tcf.csv', normalize)
  str_print = 'The velocity-velocity auto correlation function is written in %s' %(tcf_file)
  print (data_op.str_wrap(str_print, 80), flush=True)

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

  return wave_num, intensity, intensity_fit

def power_spectrum_mode(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, end_step, max_frame_corr, \
                        cluster_group_id, lower_wave, upper_wave, increment_wave, traj_coord_file, traj_vel_file, a_vec, b_vec, c_vec, normalize, work_dir):

  #Reference literature: Chem. Phys. 1986, 106, 205-212.

  '''
  power_spectrum_mode: calculate mode power_spectrum of the system

  Args :
    atoms_num: int
      atoms_num is the number of atoms in the system.
    pre_base_block: int
      pre_base_block is the number of lines before structure in a structure block.
    end_base_block: int
      end_base_block is the number of lines after structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    each: int
      each is printing frequency of md.
    start_frame_id: int
      start_frame_id is the starting frame id in the trajectory file.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    max_frame_corr: int
      max_frame_corr is the max number of correlation frames.
    cluster_group_id: 1-d int list
      cluster_group_id is the id of first atoms in the molecules in the group.
    lower_wave: float
      lower_wave is the starting wave value. The unit is cm^-1.
    upper_wave: float
      upper_wave is the ending wave value. The unit is cm^-1.
    increment_wave: float
      increment_wave is the increment of wave. The unit is cm^-1.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    traj_vel_file: string
      traj_vel_file is the name of velocity trajectory file.
    a_vec: 1-d float list, dim = 3
      a_vec is the cell vector a.
      Example : [12.42, 0.0, 0.0]
    b_vec: 1-d float list, dim = 3
      b_vec is the cell vector b.
      Example : [0.0, 12.42, 0.0]
    c_vec: 1-d float list, dim = 3
      c_vec is the cell vector c.
      Example: [0.0, 0.0, 12.42]
    normalize: int
      normalize is whether to use normalize. There are two choices: 0 and 1.
      0 means not using normalize, 1 means using normalize.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns :
    wave_num: 1-d float list
      wave_num is the list of wave number.
    Q1_intensity: 1-d float array
      intensity is the list of Q1 mode intensity for each wave number.
    Q1_intensity_fit: 1-d float array
      intensity_fit is the list of Q1 mode fitted intensity for each wave number.
    Q2_intensity: 1-d float array
      intensity is the list of Q2 mode intensity for each wave number.
    Q2_intensity_fit: 1-d float array
      intensity_fit is the list of Q2 mode fitted intensity for each wave number.
    Q3_intensity: 1-d float array
      intensity is the list of Q3 mode intensity for each wave number.
    Q3_intensity_fit: 1-d float array
      intensity_fit is the list of Q3 mode fitted intensity for each wave number.
  '''

  print ('Calculate velocity-velocity auto correlation function at first', flush=True)
  Q1_vel_tcf, Q2_vel_tcf, Q3_vel_tcf, tcf_q1_file, tcf_q2_file, tcf_q3_file = \
  time_correlation.time_corr_mode_func(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, end_step, \
                                       max_frame_corr, cluster_group_id, traj_coord_file, traj_vel_file, a_vec, b_vec, c_vec, work_dir, normalize)

  str_print = 'The velocity-velocity auto correlation function of mode 1 (symmetric strech) is written in %s' %(tcf_q1_file)
  print (data_op.str_wrap(str_print, 80), flush=True)
  str_print = 'The velocity-velocity auto correlation function of mode 2 (bending) is written in %s' %(tcf_q2_file)
  print (data_op.str_wrap(str_print, 80), flush=True)
  str_print = 'The velocity-velocity auto correlation function of mode 3 (asymmetric strech) is written in %s' %(tcf_q3_file)
  print (data_op.str_wrap(str_print, 80), flush=True)

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

  Q2_intensity = \
  statistic_mod.statistic.fourier_transform(Q2_vel_tcf, time, wave_num, time_step)
  Q2_intensity_fit = numeric.savitzky_golay(Q2_intensity,201,3)

  Q3_intensity = \
  statistic_mod.statistic.fourier_transform(Q3_vel_tcf, time, wave_num, time_step)
  Q3_intensity_fit = numeric.savitzky_golay(Q3_intensity,201,3)

  return wave_num, Q1_intensity, Q1_intensity_fit, Q2_intensity, Q2_intensity_fit, Q3_intensity, Q3_intensity_fit

def power_spectrum_run(spectrum_param, work_dir):

  '''
  power_spectrum_run: the kernel function to run power_spectrum function.

  Args:
    spectrum_param: dictionary
      spectrum_param contains keywords used in power_spectrum functions.
    work_dir: string
      work_dir is working directory of CP2K_kit.
  Returns:
    none
  '''

  spectrum_param = check_analyze.check_spectrum_inp(spectrum_param)

  spec_type = spectrum_param['type']

  if ( spec_type == 'water_mode' or spec_type == 'hydration_mode' ):
    traj_coord_file = spectrum_param['traj_coord_file']

  traj_vel_file = spectrum_param['traj_vel_file']
  init_step = int(spectrum_param['init_step'])
  end_step = int(spectrum_param['end_step'])
  max_frame_corr = int(spectrum_param['max_frame_corr'])
  start_wave = float(spectrum_param['start_wave'])
  end_wave = float(spectrum_param['end_wave'])
  increment_wave = float(spectrum_param['increment_wave'])
  normalize = int(spectrum_param['normalize'])

  if ( spec_type == 'general' ):
    atom_id = spectrum_param['atom_id']

    atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_vel_file, 'vel')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    print ('POWER_SPECTRUM'.center(80, '*'), flush=True)
    print ('Analyze the power spectrum for the choosed system', flush=True)
    wave_num, intensity, intensity_fit = power_spectrum(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, \
                                         end_step, max_frame_corr, start_wave, end_wave, increment_wave, atom_id, traj_vel_file, normalize, work_dir)

    freq_int_file = ''.join((work_dir, '/freq_intensity.csv'))
    with open(freq_int_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      if ( normalize == 0 ):
        writer.writerow(['wave(cm^-1)', 'intensity(cm^2/s)', 'intensity_fit(cm^2/s)'])
      elif ( normalize == 1):
        writer.writerow(['wave(cm^-1)', 'intensity', 'intensity_fit'])
      for i in range(len(wave_num)):
        writer.writerow([wave_num[i], intensity[i], intensity_fit[i]])
    str_tmp = 'The power spectrum is written in %s' %(freq_int_file)
    print (data_op.str_wrap(str_tmp, 80), flush=True)

  else:
    atoms_num_p, pre_base_block_p, end_base_block_p, pre_base_p, frames_num_p, each_p, start_frame_id_p, end_frame_id_p, time_step_p = \
    traj_info.get_traj_info(traj_coord_file, 'coord_xyz')
    atoms_num_v, pre_base_block_v, end_base_block_v, pre_base_v, frames_num_v, each_v, start_frame_id_v, end_frame_id_v, time_step_v = \
    traj_info.get_traj_info(traj_vel_file, 'vel')

    a_vec = spectrum_param['box']['A']
    b_vec = spectrum_param['box']['B']
    c_vec = spectrum_param['box']['C']

    if ( spec_type == 'water_mode' ):
      atom_id = spectrum_param['atom_id']
      if (start_frame_id_v > start_frame_id_p):
        atoms_num, pre_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_vel_file, 'vel')
      else:
        atoms_num, pre_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_coord_file, 'coord_xyz')

      log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

      if ( len(atom_id) != atoms_num ):
        choose_coord_file = traj_tools.choose_str(atoms_num, pre_base, pre_base_block, end_base_block, each, start_frame_id, \
                            end_frame_id, start_frame_id, traj_coord_file, atom_id, work_dir, 'choose_coord.xyz')
        choose_vel_file = traj_tools.choose_str(atoms_num, pre_base, pre_base_block, end_base_block, each, start_frame_id, \
                            end_frame_id, start_frame_id, traj_vel_file, atom_id, work_dir, 'choose_coord.xyz')
        atoms_num = len(atom_id)
        atom_id = list(range(1,atoms_num+1,1))
        coord_order_file, order_list = \
        geometry.order_struct(atoms_num, frames_num, pre_base_block, end_base_block, pre_base, [['O','H','H']], \
                              [atom_id], choose_coord_file, a_vec, b_vec, c_vec, work_dir, 'coord_order.xyz')
        vel_order_file = traj_tools.order_traj_file(atoms_num, frames_num, each, start_frame_id, \
                         choosed_vel_file, 'vel', order_list, work_dir, 'vel_order.xyz')
      else:
        coord_order_file, order_list = \
        geometry.order_struct(atoms_num, frames_num, pre_base_block, end_base_block, pre_base, [['O','H','H']], \
                              [atom_id], traj_coord_file, a_vec, b_vec, c_vec, work_dir, 'coord_order.xyz')
        vel_order_file = traj_tools.order_traj_file(atoms_num, frames_num, each, start_frame_id, \
                         traj_vel_file, 'vel',  order_list, work_dir, 'vel_order.xyz')

      traj_coord_file = coord_order_file
      traj_vel_file = vel_order_file

      if (start_frame_id_v > start_frame_id_p):
        atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, \
        time_step, group_atom_1_id, group_atoms_mass = \
        traj_info.get_traj_info(traj_vel_file, 'vel', [['O','H','H']], [atom_id], True)
      else:
        atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, \
        time_step, group_atom_1_id, group_atoms_mass = \
        traj_info.get_traj_info(traj_coord_file, 'coord_xyz', [['O','H','H']], [atom_id], True)

      cluster_id = []
      for i in range(int(end_step)-int(init_step)+1):
        cluster_frame_id = []
        for j in range(len(group_atom_1_id[0])):
          cluster_frame_id.append([group_atom_1_id[0][j], group_atom_1_id[0][j]+1, group_atom_1_id[0][j]+2])
        cluster_id.append(cluster_frame_id)

      print ('POWER_SPECTRUM'.center(80, '*'), flush=True)
      print ('Analyze power spectrum of water with three modes of water', flush=True)
      wave_num, q1_int, q1_int_fit, q2_int, q2_int_fit, q3_int, q3_int_fit = \
      power_spectrum_mode(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, end_step, max_frame_corr, \
                          cluster_id, start_wave, end_wave, increment_wave, traj_coord_file, traj_vel_file, a_vec, b_vec, c_vec, normalize, work_dir)

    elif ( spec_type == 'hydration_mode' ):
      if (start_frame_id_v > start_frame_id_p):
        atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_id, time_step = \
        traj_info.get_traj_info(traj_vel_file, 'vel')
      else:
        atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_id, time_step = \
        traj_info.get_traj_info(traj_coord_file, 'coord_xyz')

      log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_id, time_step)

      hyd_shell_dist = spectrum_param['hyd_shell_dist']
      dist_conv = spectrum_param['dist_conv']
      atom_type_pair = spectrum_param['atom_type_pair']
      atom_1 = atom_type_pair[0]
      atom_2 = atom_type_pair[1]

      first_shell_id, dist = geometry.first_shell(atoms_num, pre_base_block, end_base_block, pre_base, start_frame_id, frames_num, each, \
                             init_step, end_step, atom_1, atom_2, a_vec, b_vec, c_vec, traj_coord_file, hyd_shell_dist, dist_conv, work_dir)

      cluster_id = []
      frame_num_stat = int((end_step-init_step)/each+1)
      for i in range(frame_num_stat):
        frame_id = init_step+i*each
        angle_order = geometry.order_angle(first_shell_id[i][0][0], first_shell_id[i][0][1:len(first_shell_id[i][0])], \
                                           frame_id, each, traj_coord_file)
        angle_order_cyclic = angle_order+angle_order
        cluster_frame_id = []
        for j in range(len(angle_order)):
          group_id = []
          group_id.append(first_shell_id[i][0][0])
          group_id.append(angle_order_cyclic[j])
          group_id.append(angle_order_cyclic[j+1])
          cluster_frame_id.append(group_id)
        cluster_id.append(cluster_frame_id)

      print ('POWER_SPECTRUM'.center(80, '*'), flush=True)
      print ('Analyze power spectrum for metal ion and its first shell water', flush=True)

      wave_num, q1_int, q1_int_fit, q2_int, q2_int_fit, q3_int, q3_int_fit = \
      power_spectrum_mode(atoms_num, pre_base_block, end_base_block, pre_base, each, start_frame_id, time_step, init_step, end_step, max_frame_corr, \
                          cluster_id, start_wave, end_wave, increment_wave, traj_coord_file, traj_vel_file, a_vec, b_vec, c_vec, normalize, work_dir)

    freq_int_q1_file = ''.join((work_dir, '/freq_intensity_q1.csv'))
    with open(freq_int_q1_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      if ( normalize == 0 ):
        writer.writerow(['wave(cm^-1)', 'intensity(cm^2/s)', 'intensity_fit(cm^2/s)'])
      elif ( normalize == 1 ):
        writer.writerow(['wave(cm^-1)', 'intensity', 'intensity_fit'])
      for i in range(len(wave_num)):
        writer.writerow([wave_num[i], q1_int[i], q1_int_fit[i]])

    freq_int_q2_file = ''.join((work_dir, '/freq_intensity_q2.csv'))
    with open(freq_int_q2_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      if ( normalize == 0 ):
        writer.writerow(['wave(cm^-1)', 'intensity(cm^2/s)', 'intensity_fit(cm^2/s)'])
      elif ( normalize == 1 ):
        writer.writerow(['wave(cm^-1)', 'intensity', 'intensity_fit'])
      for i in range(len(wave_num)):
        writer.writerow([wave_num[i], q2_int[i], q2_int_fit[i]])

    freq_int_q3_file = ''.join((work_dir, '/freq_intensity_q3.csv'))
    with open(freq_int_q3_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      if ( normalize == 0 ):
        writer.writerow(['wave(cm^-1)', 'intensity(cm^2/s)', 'intensity_fit(cm^2/s)'])
      elif ( normalize == 1 ):
        writer.writerow(['wave(cm^-1)', 'intensity', 'intensity_fit'])
      for i in range(len(wave_num)):
        writer.writerow([wave_num[i], q3_int[i], q3_int_fit[i]])

    str_print = 'The power spectrum of mode 1 (symmetric strech) is written in %s' %(freq_int_q1_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    str_print = 'The power spectrum of mode 2 (bending) is written in %s' %(freq_int_q2_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    str_print = 'The power spectrum of mode 3 (asymmetric) is written in %s' %(freq_int_q3_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

