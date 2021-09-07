#!/usr/bin/env python

import os
import csv
import linecache
import numpy as np
from CP2K_kit.tools import atom
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.lib import dynamic_mod
from CP2K_kit.analyze import check_analyze

def diffusion_msd(atoms_num, base, pre_base, each, start_frame_id, time_step, init_step, end_step, \
                  max_frame_corr, atom_id, traj_coord_file, remove_com, work_dir, file_name):

  '''
  diffusion_msd: calculate diffusion coefficient for choosed atoms via mean square displacement.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of trajectory.
    each: int
      each is printing frequency of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in trajectory file.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    max_frame_corr: int
      max_frame_corr is the max number of correlation frames.
    atom_id: 1-d int list
      atom_id is the id of atoms.
    traj_coord_file: string
      file_name is the name of coordination trajectory file.
    remove_com: bool
      remove_com is whether we need to remove center of mass.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    msd_file: string
      msd_file is the file_name of mean square displacement.
  '''

  frames_num_stat = int((end_step-init_step)/each+1)
  coord = np.asfortranarray(np.zeros((frames_num_stat, len(atom_id), 3)), dtype='float32')

  #Dump coordinate
  for i in range(frames_num_stat):
    for j in range(len(atom_id)):
      line_ij = linecache.getline(traj_coord_file, (atoms_num+base)*(int((init_step-start_frame_id)/each)+i)+pre_base+base+atom_id[j])
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      coord[i,j,0] = float(line_ij_split[1])
      coord[i,j,1] = float(line_ij_split[2])
      coord[i,j,2] = float(line_ij_split[3])

  if remove_com:
    #Dump atom mass
    atom_mass = []
    for i in range(len(atom_id)):
      line_i = linecache.getline(traj_coord_file, atom_id[i]+pre_base+base)
      line_i_split = data_op.split_str(line_i, ' ')
      atom_mass.append(atom.get_atom_mass(line_i_split[0])[1])
    atom_mass_array = np.asfortranarray(atom_mass, dtype='float32')
    coord = dynamic_mod.dynamic.remove_coord_com(coord,atom_mass_array)

  linecache.clearcache()

  einstein_sum = dynamic_mod.dynamic.diffusion_einstein_sum(coord, max_frame_corr)

  msd_file = ''.join((work_dir, '/', file_name))
  with open(msd_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time(fs)', 'msd(Angstrom^2)'])
    for i in range(len(einstein_sum)):
      writer.writerow([i*time_step*each, einstein_sum[i]])

  return msd_file

def diffusion_tcf(atoms_num, base, pre_base, each, start_frame_id, time_step, \
                  init_step, end_step, max_frame_corr, atom_id, traj_vel_file):

  '''
  diffusion_tcf: calculate diffusion coefficient for choosed atoms via velocity time correlation function.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in trajectory file.
    base: int
      base is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of trajectory.
    each: int
      each is printing frequency of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in trajectory file.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the ending step frame id.
    max_frame_corr: int
      max_frame_corr is the max number of correlation frames.
    atom_id: 1-d int list
      atom_id is the id of atoms.
    traj_vel_file: string
      traj_vel_file is the name of velocity trajectory file.
  Returns:
    diff_coeff: float
      diff_coeff is the diffusion coefficient.
  '''

  #Do we need to substract com velocity?

  frames_num_stat = int((end_step-init_step)/each+1)
  vel = np.asfortranarray(np.zeros((frames_num_stat, len(atom_id), 3)), dtype='float32')
  #Dump velocity
  for i in range(frames_num_stat):
    for j in range(len(atom_id)):
      line_ij = linecache.getline(traj_vel_file, (atoms_num+base)*(int((init_step-start_frame_id)/each)+i)+base+atom_id[j])
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      vel[i,j,0] = float(line_ij_split[1])
      vel[i,j,1] = float(line_ij_split[2])
      vel[i,j,2] = float(line_ij_split[3])

  linecache.clearcache()

  #Here we use non-normalized velocity time correlation function.
  normalize = 0

  vel_tcf = dynamic_mod.dynamic.time_correlation(vel, max_frame_corr, normalize)

  sum_value = 0.0
  for i in range(len(vel_tcf)):
    sum_value = sum_value+vel_tcf[i]*time_step*each*1.0E-15
  diff_coeff = sum_value/3.0

  return diff_coeff

def diffusion_run(diffusion_param, work_dir):

  '''
  diffusion_run: kernel function to run diffusion. It will call diffusion function.

  Args:
    diffusion_param: dictionary
      diffusion_param contains keywords used in diffusion function.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  diffusion_param = check_analyze.check_diffusion_inp(diffusion_param)

  method = diffusion_param['method']
  atom_id = diffusion_param['atom_id']
  init_step = diffusion_param['init_step']
  end_step = diffusion_param['end_step']
  max_frame_corr = diffusion_param['max_frame_corr']

  if ( method == 'einstein_sum' ):
    traj_coord_file = diffusion_param['traj_coord_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_coord_file, 'coord')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    print ('DIFFUSION'.center(80, '*'), flush=True)
    print ('Diffusion coefficient is calculated by %s' %(method), flush=True)

    remove_com = diffusion_param['remove_com']
    msd_file = diffusion_msd(atoms_num, base, pre_base, each, start_frame_id, time_step, init_step, \
               end_step, max_frame_corr, atom_id, traj_coord_file, remove_com, work_dir, 'msd.csv')

    str_tmp = 'The mean square displacement file is written in %s' %(msd_file)
    print (data_op.str_wrap(str_tmp, 80), flush=True)
    str_tmp = 'To get the diffusion coefficient, user should get the slope of msd vs time, and\
    diffusion coefficient is slope/6.0*0.1 cm^2/s'
    print (data_op.str_wrap(str_tmp, 80), flush=True)

  elif ( method == 'green_kubo' ):
    traj_vel_file = diffusion_param['traj_vel_file']
    atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
    traj_info.get_traj_info(traj_vel_file, 'vel')

    log_info.log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step)

    print ('DIFFUSION'.center(80, '*'), flush=True)
    print ('Diffusion coefficient is calculated by %s' %(method), flush=True)

    diff_coeff = diffusion_tcf(atoms_num, base, pre_base, each, start_frame_id, time_step, init_step, \
                 end_step, max_frame_corr, atom_id, traj_vel_file)

    print ("The diffusion coefficient calculated by vel_tcf is %f cm^2/s" %diff_coeff, flush=True)

