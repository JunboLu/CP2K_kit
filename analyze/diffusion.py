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

def diffusion(atoms_num, base, pre_base, each, file_start, time_step, start, end, \
              time_num, atom_id, file_name, method, remove_com, work_dir):

  '''
  diffusion : calculate diffusion coefficient for choosed atoms.

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
    method : string
      method is the method to calculate diffusion coefficient. Two method choices: einstein and green-Kubo.
      The einstein mthod is better!
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  '''
  For einstein method, we should be careful for unit:
  The unit of MSD is Angstrom^2
  The unit of time is fs
  The diffusion coeff is slope/6.0*0.1
  '''

  frames_num_stat = int((end-start)/each+1)
  if (method == 'einstein_sum'):
    coord = np.asfortranarray(np.zeros((frames_num_stat, len(atom_id), 3)), dtype='float32')

    #Dump coordinate
    for i in range(frames_num_stat):
      for j in range(len(atom_id)):
        line_ij = linecache.getline(file_name, (atoms_num+base)*(int((start-file_start)/each)+i)+pre_base+base+atom_id[j])
        line_ij_split = data_op.str_split(line_ij, ' ')
        coord[i,j,0] = float(line_ij_split[1])
        coord[i,j,1] = float(line_ij_split[2])
        coord[i,j,2] = float(line_ij_split[3].strip('\n'))

    if remove_com:
      #Dump atom mass
      atom_mass = []
      for i in range(len(atom_id)):
        line_i = linecache.getline(file_name, atom_id[i]+pre_base+base)
        line_i_split = data_op.str_split(line_i, ' ')
        atom_mass.append(atom.get_atom_mass(line_i_split[0])[1])
      atom_mass_array = np.asfortranarray(atom_mass, dtype='float32')
      coord = dynamic_mod.dynamic.remove_coord_com(coord,atom_mass_array)

    einstein_sum = dynamic_mod.dynamic.diffusion_einstein_sum(coord, time_num)

    msd_file = ''.join((work_dir, '/msd.csv'))
    with open(msd_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time(ps)', 'msd(Angstrom^2)'])
      for i in range(len(einstein_sum)):
        writer.writerow([i*time_step*each*0.001, einstein_sum[i]])

  if (method == 'green_kubo'):
    vel = np.asfortranarray(np.zeros((frames_num_stat, len(atom_id), 3)), dtype='float32')
    #Dump velocity
    for i in range(frames_num_stat):
      for j in range(len(atom_id)):
        line_ij = linecache.getline(file_name, (atoms_num+base)*(int((start-file_start)/each)+i)+base+atom_id[j])
        line_ij_split = data_op.str_split(line_ij, ' ')
        vel[i,j,0] = float(line_ij_split[1])
        vel[i,j,1] = float(line_ij_split[2])
        vel[i,j,2] = float(line_ij_split[3].strip('\n'))
    #Here we use non-normalized velocity time correlation function.
    normalize = 0

    vel_tcf = dynamic_mod.dynamic.time_correlation(vel, time_num, normalize)

    sum_value = 0.0
    for i in range(len(vel_tcf)):
      sum_value = sum_value+vel_tcf[i]*time_step*each*1.0E-15
    Diff_coeff = sum_value/3.0
    diff_file_name = ''.join((work_dir, '/diffusion'))
    diff_file = open(diff_file_name, 'w')
    diff_file.write("The diffusion coefficient by vel_tcf is %f cm^2/s" %Diff_coeff)
    diff_file.close()

def diffusion_run(diffusion_param, work_dir):

  '''
  diffusion_run : kernel function to run diffusion. It will call diffusion function.

  Args :
    diffusion_param : dictionary
      diffusion_param contains keywords used in diffusion function.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  if ( 'traj_file' in diffusion_param.keys() ):
    traj_file = diffusion_param['traj_file']
    if ( os.path.exists(traj_file) ):
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(traj_file)
    else:
      log_info.log_error('%s file does not exist' %(traj_file))
  else:
    log_info.log_error('No trajectory file found, please set analyze/diffusion/traj_file')
    exit()

  if ( 'method' in diffusion_param.keys() ):
    method = diffusion_param['method']
  else:
    method = 'einstein_sum'

  if ( 'atom_id' in diffusion_param.keys() ):
    atom_id = data_op.get_id_list(diffusion_param['atom_id'])
  else:
    log_info.log_error('No atom id found, please set analyze/diffusion/atom_id')
    exit()

  if ( 'init_step' in diffusion_param.keys() ):
    init_step = int(diffusion_param['init_step'])
  else:
    init_step = start_frame_id

  if ( 'end_step' in diffusion_param.keys() ):
    end_step = int(diffusion_param['end_step'])
  else:
    end_step = start_frame_id + each

  if ( 'max_frame_corr' in diffusion_param.keys() ):
    max_frame_corr = int(diffusion_param['max_frame_corr'])
  else:
    max_frame_corr = int(frame_num/2)

  if ( 'remove_com' in diffusion_param.keys() ):
    remove_com = data_op.str_to_bool(diffusion_param['remove_com'])
  else:
    remove_com = False

  diffusion(atoms_num, base, pre_base, each, start_frame_id, time_step, init_step, \
            end_step, max_frame_corr, atom_id, traj_file, method, remove_com, work_dir)
