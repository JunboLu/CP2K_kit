#! /usr/env/bin python

import os
import copy
import math
import linecache
import numpy as np
from CP2K_kit.tools import call
from CP2K_kit.tools import numeric
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import log_info

hartree_to_ev = 2.72113838565563E+01
ang_to_bohr = 1.0/5.29177208590000E-01

def load_data_from_sepfile(file_dir, file_prefix, proj_name, tot_atoms_type_dic):

  '''
  load_data_from_sepfile: load training data from separate files.

  Args:
    file_dir: string
      file_dir is the directory of initial training data.
    file_prefix: string
      file_prefix is the prefix of file name.
    proj_name: string
      proj_name is the cp2k project name.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  cmd = "ls | grep %s" %(file_prefix)
  task_num = len(call.call_returns_shell(file_dir, cmd))

  if ( task_num != 0 ):
    cmd = "mkdir %s" % ('data')
    call.call_simple_shell(file_dir, cmd)
    save_dir = ''.join((file_dir, '/data'))
    type_file = open(''.join((save_dir, '/type.raw')), 'w')
    box_file = open(''.join((save_dir, '/box.raw')), 'w')
    frc_file = open(''.join((save_dir, '/force.raw')), 'w')
    coord_file = open(''.join((save_dir, '/coord.raw')), 'w')
    energy_file = open(''.join((save_dir, '/energy.raw')), 'w')
    virial_file = open(''.join((save_dir, '/virial.raw')), 'w')

    atoms = []

    coord_file_0 = ''.join((file_dir, '/', file_prefix, str(0), '/coord'))
    atoms_num = len(open(coord_file_0).readlines())

    #Dump atoms type information
    for i in range(atoms_num):
      line_i = linecache.getline(coord_file_0, i+1)
      line_i_split = data_op.split_str(line_i, ' ')
      atoms.append(line_i_split[0])

    linecache.clearcache()

    for i in range(atoms_num):
      atom_type_index = tot_atoms_type_dic[atoms[i]]
      type_file.write(''.join((str(atom_type_index), '\n')))

    type_file.close()

    for i in range(task_num):

      task_dir_i = ''.join((file_dir, '/', file_prefix, str(i)))
      box_file_i = ''.join((task_dir_i, '/box'))
      coord_file_i = ''.join((task_dir_i, '/coord'))
      out_file_i = ''.join((task_dir_i, '/cp2k.out'))
      frc_file_i = ''.join((task_dir_i, '/', proj_name, '-1_0.xyz'))
      stress_file_i = ''.join((task_dir_i, '/', proj_name, '-1_0.stress_tensor'))

      box_exist = os.path.exists(box_file_i)
      coord_exist = os.path.exists(coord_file_i)
      out_exist = os.path.exists(out_file_i)
      frc_exist = os.path.exists(frc_file_i)

      #Dump box information
      if ( box_exist and coord_exist and out_exist and frc_exist ):
        grep_scf = 'grep %s %s' %("'SCF run NOT converged'", out_file_i)
        grep_scf_info = call.call_returns_shell(save_dir, grep_scf)
        if ( len(grep_scf_info) == 0 ):
          vec = []
          for j in range(3):
            line_ij = linecache.getline(box_file_i, j+1)
            line_ij_split = data_op.split_str(line_ij, ' ', '\n')
            #vol will be used in stress calculation.
            vec.append([float(line_ij_split[1]), \
                        float(line_ij_split[2]), \
                        float(line_ij_split[3])])
            if ( j== 0 ):
              frame_str = ' '.join((line_ij_split[1], line_ij_split[2], line_ij_split[3]))
            else:
              frame_str = ' '.join((frame_str, line_ij_split[1], line_ij_split[2], line_ij_split[3]))

          linecache.clearcache()

          vol = np.linalg.det(vec)
          frame_str = ''.join((frame_str, '\n'))
          box_file.write(frame_str)

          #Dump energy information
          cmd = "grep %s %s" % ("'ENERGY| Total FORCE_EVAL'", 'cp2k.out')
          energy_parse = call.call_returns_shell(task_dir_i, cmd)
          energy = float(energy_parse[0].split(':')[1].strip('\n'))
          energy_ev = energy*hartree_to_ev
          energy_file.write(''.join((str(energy_ev), '\n')))

          #Dump coord information
          for j in range(atoms_num):
            line_ij = linecache.getline(coord_file_i, j+1)
            line_ij_split = data_op.split_str(line_ij, ' ', '\n')
            if (j==0):
              frame_str = ' '.join((line_ij_split[1], line_ij_split[2], line_ij_split[3]))
            else:
              frame_str = ' '.join((frame_str, line_ij_split[1], line_ij_split[2], line_ij_split[3]))

          linecache.clearcache()

          frame_str = ''.join((frame_str, '\n'))
          coord_file.write(frame_str)

          #Dump force information
          for j in range(atoms_num):
            line_ij = linecache.getline(frc_file_i, j+4+1)
            line_ij_split = data_op.split_str(line_ij, ' ', '\n')
            f1 = float(line_ij_split[3])*hartree_to_ev*ang_to_bohr
            f2 = float(line_ij_split[4])*hartree_to_ev*ang_to_bohr
            f3 = float(line_ij_split[5])*hartree_to_ev*ang_to_bohr
            if (j==0):
              frame_str = ' '.join((str(f1), str(f2), str(f3)))
            else:
              frame_str = ' '.join((frame_str, str(f1), str(f2), str(f3)))

          linecache.clearcache()

          frame_str = ''.join((frame_str, '\n'))
          frc_file.write(frame_str)

          #Dump virial information
          #stress in 'xx-1.stress_tensor' is GPa.
          if os.path.exists(stress_file_i):
            for j in range(3):
              line_ij = linecache.getline(stress_file_i, j+5)
              line_ij_split = data_op.split_str(line_ij, ' ', '\n')
              stress_1 = float(line_ij_split[1])*vol/160.21766208
              stress_2 = float(line_ij_split[2])*vol/160.21766208
              stress_3 = float(line_ij_split[3])*vol/160.21766208
              stress_str_1 = numeric.get_as_num_string(stress_1)
              stress_str_2 = numeric.get_as_num_string(stress_2)
              stress_str_3 = numeric.get_as_num_string(stress_3)
              if ( j==0 ):
                frame_str = ' '.join((stress_str_1, stress_str_2, stress_str_3))
              else:
                frame_str = ' '.join((frame_str, stress_str_1, stress_str_2, stress_str_3))

            linecache.clearcache()

            frame_str = ''.join((frame_str, '\n'))
            virial_file.write(frame_str)

    box_file.close()
    frc_file.close()
    coord_file.close()
    energy_file.close()
    virial_file.close()

    if ( os.path.getsize(''.join((save_dir, '/virial.raw'))) == 0 ):
      cmd = "rm %s" % (''.join((save_dir, '/virial.raw')))
      call.call_simple_shell(save_dir, cmd)

def load_data_from_dir(traj_coord_file_name, traj_frc_file_name, traj_cell_file_name, traj_stress_file_name, \
                       train_stress, work_dir, save_dir, start, end, choosed_num, tot_atoms_type_dic):

  '''
  load_data_from_dir: load initial training data from cp2k calculation directory.

  Args:
    traj_coord_file_name: string
      traj_coord_file_name is the name of coordination trajectory file.
    traj_frc_file_name: string
      traj_frc_file_name is the name of force trajectory file.
    traj_cell_file_name: string
      traj_cell_file_name is the name of cell trajectory file.
    traj_stress_file_name: string
      traj_stress_file_name is the name of stress trajectory file.
    train_stress: bool
      train_stress is whether we need to train stress
    work_dir: string
      work_dir is workding directory.
    save_dir: string
      save_dir is the directory where data will be saved.
    start: int
      start is the starting id of trajectory.
    end: int
      end is the endding id of trajectory.
    choosed_num: int
      choosed_num is the number of choosed frames.
    tot_atoms_type_dic: dictionary
      tot_atoms_type_dic is the atoms type dictionary.
  Returns:
    none
  '''

  cmd = "mkdir %s" % (save_dir)
  call.call_simple_shell(work_dir, cmd)

  atoms_num, base, pre_base, frames_num, each, start_id, end_id, time_step = \
  traj_info.get_traj_info(traj_coord_file_name, 'coord')

  if ( start < start_id ):
    log_info.log_error('Input error: start is less than the starting frame id in trajectory, please check or reset deepff/deepmd/training/system/start')
    exit()

  if ( end > end_id ):
    log_info.log_error('Input error: end is larger than the endding frame id in trajectory, please check or reset deepff/deepmd/training/system/end')
    exit()

  if ( start >= end ):
    log_info.log_error('Input error: start_frame is larger than end_frame, please check!')
    exit()
  else:
    total_index = data_op.gen_list(start, end, each)
    max_choosed_num = int((end-start)/each)+1
    if ( choosed_num > max_choosed_num ):
      log_info.log_error('Input error: choosed_frame_num is larger than the number of frames in trajectory, please check!')
      exit()
    else:
      total_index_array = np.array(total_index)
      if ( choosed_num < max_choosed_num ):
        np.random.shuffle(total_index_array)
      choosed_index = list(total_index_array[0:choosed_num])

  #Dump box information
  box_file = open(''.join((save_dir, '/box.raw')),'w')
  #vol is volume, it will be used in virial. The unit of vol is A^3.
  vol = []
  for i in range(len(choosed_index)):
    line_i = linecache.getline(traj_cell_file_name, int((choosed_index[i]-start_id)/each)+2)
    line_i_split = data_op.split_str(line_i, ' ', '\n')
    vol.append(float(line_i_split[len(line_i_split)-1]))
    for j in range(9):
      if ( j == 0 ):
        frame_str = ''.join((line_i_split[j+2]))
      else:
        frame_str = ' '.join((frame_str, line_i_split[j+2]))
    frame_str = ''.join((frame_str, '\n'))
    box_file.write(frame_str)

  linecache.clearcache()

  box_file.close()

  #Dump energy information
  energy_file = open(''.join((save_dir, '/energy.raw')), 'w')
  for i in range(len(choosed_index)):
    line_i = linecache.getline(traj_coord_file_name, int((choosed_index[i]-start_id)/each)*(atoms_num+2) + 2)
    line_i_split = data_op.split_str(line_i, ' ', '\n')
    energy = float(line_i_split[len(line_i_split)-1])*hartree_to_ev
    energy_file.write(''.join((str(energy),'\n')))

  linecache.clearcache()

  #Dump element information
  type_file = open(''.join((save_dir, '/type.raw')), 'w')
  atoms = []
  for i in range(atoms_num):
    line_i = linecache.getline(traj_coord_file_name, 2+i+1)
    line_i_split = data_op.split_str(line_i, ' ')
    atoms.append(line_i_split[0])

  linecache.clearcache()

  for i in range(atoms_num):
    atom_type_index = tot_atoms_type_dic[atoms[i]]
    type_file.write(''.join((str(atom_type_index), '\n')))

  #Dump coordination information
  coord_file = open(''.join((save_dir, '/coord.raw')), 'w')
  for i in range(len(choosed_index)):
    frame_str = ''
    for j in range(atoms_num):
      line_ij = linecache.getline(traj_coord_file_name, int((choosed_index[i]-start_id)/each)*(atoms_num+2)+2+j+1)
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      if (j==0):
        frame_str = ' '.join((line_ij_split[1], line_ij_split[2], line_ij_split[3]))
      else:
        frame_str = ' '.join((frame_str, line_ij_split[1], line_ij_split[2], line_ij_split[3]))

    frame_str = ''.join((frame_str, '\n'))
    coord_file.write(frame_str)

  linecache.clearcache()
  coord_file.close()

  #Dump force information
  frc_file = open(''.join((save_dir, '/force.raw')), 'w')
  for i in range(len(choosed_index)):
    frame_str = ''
    for j in range(atoms_num):
      line_ij = linecache.getline(traj_frc_file_name, int((choosed_index[i]-start_id)/each)*(atoms_num+2)+2+j+1)
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      f1 = float(line_ij_split[1])*hartree_to_ev*ang_to_bohr
      f2 = float(line_ij_split[2])*hartree_to_ev*ang_to_bohr
      f3 = float(line_ij_split[3])*hartree_to_ev*ang_to_bohr
      if ( j == 0 ):
        frame_str = ' '.join((str(f1), str(f2), str(f3)))
      else:
        frame_str = ' '.join((frame_str, str(f1), str(f2), str(f3)))

    frame_str = ''.join((frame_str, '\n'))
    frc_file.write(frame_str)

  linecache.clearcache()

  frc_file.close()

  #Dump virial information
  if ( train_stress and traj_stress_file_name != 'none' ):
    virial_file = open(''.join((save_dir, '/virial.raw')), 'w')
    for i in range(len(choosed_index)):
      frame_str = ''
      #There are 9 elements: xx, xy, xz, yx, yy, yz, zx, zy, zz.
      line_i = linecache.getline(traj_stress_file_name, int((choosed_index[i]-start_id)/each)+2)
      line_i_split = data_op.split_str(line_i, ' ', '\n')
      #The unit of stress tensor in 'xx-1.stress' file is bar.
      for j in range(9):
        stress_j = float(line_i_split[j+2])*vol[i]/(1602176.6208)
        stress_j_str = numeric.get_as_num_string(stress_j)
        if ( j == 0 ):
          frame_str = ''.join((stress_j_str))
        else:
          frame_str = ' '.join((frame_str, stress_j_str))
      frame_str = ''.join((frame_str, '\n'))
      virial_file.write(frame_str)

    linecache.clearcache()

    virial_file.close()

  total_index = data_op.gen_list(start_id, end_id, each)
  no_train_data_index = copy.deepcopy(total_index)

  for i in choosed_index:
    no_train_data_index.remove(i)

  train_data_dir = ''.join((save_dir, '/train_data'))
  cmd = "mkdir %s" %('train_data')
  call.call_simple_shell(save_dir, cmd)
  train_pos_file_name = ''.join((train_data_dir, '/train-pos-1.xyz'))
  train_frc_file_name = ''.join((train_data_dir, '/train-frc-1.xyz'))
  train_pos_file = open(train_pos_file_name, 'w')
  train_frc_file = open(train_frc_file_name, 'w')
  for i in choosed_index:
    line_1 = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+1)
    line_2 = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+2)
    train_pos_file.write(line_1)
    train_pos_file.write(line_2)
    train_frc_file.write(line_1)
    train_frc_file.write(line_2)
    for j in range(atoms_num):
      line_pos = linecache.getline(traj_coord_file_name, int((i-start_id)/each)*(atoms_num+2)+j+1+2)
      line_frc = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+j+1+2)
      train_pos_file.write(line_pos)
      train_frc_file.write(line_frc)
  train_pos_file.close()
  train_pos_file.close()

  if ( no_train_data_index != [] ):
    no_train_data_dir = ''.join((save_dir, '/no_train_data'))
    cmd = "mkdir %s" %('no_train_data')
    call.call_simple_shell(save_dir, cmd)
    no_train_pos_file_name = ''.join((no_train_data_dir, '/no_train-pos-1.xyz'))
    no_train_frc_file_name = ''.join((no_train_data_dir, '/no_train-frc-1.xyz'))
    no_train_pos_file = open(no_train_pos_file_name, 'w')
    no_train_frc_file = open(no_train_frc_file_name, 'w')
    for i in no_train_data_index:
      line_1 = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+1)
      line_2 = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+2)
      no_train_pos_file.write(line_1)
      no_train_pos_file.write(line_2)
      no_train_frc_file.write(line_1)
      no_train_frc_file.write(line_2)
      for j in range(atoms_num):
        line_pos = linecache.getline(traj_coord_file_name, int((i-start_id)/each)*(atoms_num+2)+j+1+2)
        line_frc = linecache.getline(traj_frc_file_name, int((i-start_id)/each)*(atoms_num+2)+j+1+2)
        no_train_pos_file.write(line_pos)
        no_train_frc_file.write(line_frc)
    no_train_pos_file.close()
    no_train_frc_file.close()

  linecache.clearcache()

def read_raw_data(data_dir):

  '''
  read_raw_data: read raw data from data directory

  Args :
    data_dir: string
      data_dir is the directory of initial train data.
  Returns:
    energy_array: 1-d float array, dim = num of frames
      energy_array is the 1-d array of energies.
    coord_array: 2-d float array, dim = (num of frames)*(3*(num of atoms))
      coord_array is the 2-d array of coordinates
    frc_array: 2-d float array, dim = (num of frames)*(3*(num of atoms))
      frc_array is the 2-d array of forces.
    box_array: 2-d float array, dim = (num of frames)*9
      box_array is the 2-d array of cells.
    virial_array: 2-d float array, dim = (num of frames)*9
      virial_array is the 2-d array of virials
  '''

  #Read raw data
  ene_file = ''.join((data_dir, '/energy.raw'))
  coord_file = ''.join((data_dir, '/coord.raw'))
  frc_file = ''.join((data_dir, '/force.raw'))
  box_file = ''.join((data_dir, '/box.raw'))
  virial_file = ''.join((data_dir, '/virial.raw'))

  ene_raw = open(ene_file, "rb").read()
  coord_raw = open(coord_file, "rb").read()
  frc_raw = open(frc_file, "rb").read()
  box_raw = open(box_file, "rb").read()

  frames_num = len(ene_raw.split())
  if ( frames_num == 0 ):
    log_info.log_error('Dump new cp2k data error: number of new data is zero, check the scf converge of cp2k force calculation!')
    exit()
  else:
    atoms_num = int(len(coord_raw.split())/frames_num/3)

    ene_file_exist = os.path.exists(ene_file)
    coord_file_exist = os.path.exists(coord_file)
    frc_file_exist = os.path.exists(frc_file)
    box_file_exist = os.path.exists(box_file)
    #Organize raw data
    if ( ene_file_exist and coord_file_exist and frc_file_exist and box_file_exist ):
      ene_list = []
      ene_split = ene_raw.split()
      for i in range(len(ene_split)):
        ene_list.append(float(ene_split[i].decode()))
      energy_array = np.array(ene_list)

      coord_list = []
      coord_split = coord_raw.split()
      for i in range(len(coord_split)):
        coord_list.append(float(coord_split[i].decode()))
      coord = np.array(coord_list)
      coord_array = coord.reshape(frames_num, atoms_num*3)

      frc_list = []
      frc_split = frc_raw.split()
      for i in range(len(frc_split)):
        frc_list.append(float(frc_split[i].decode()))
      frc = np.array(frc_list)
      frc_array = frc.reshape(frames_num, atoms_num*3)

      box_list = []
      box_split = box_raw.split()
      for i in range(len(box_split)):
        box_list.append(float(box_split[i].decode()))
      box = np.array(box_list)
      box_array = box.reshape(frames_num, 9)

      #virial is not necessary
      if ( os.path.exists(virial_file) ):
        virial_raw = open(virial_file, "rb").read()
        virial_list = []
        virial_split = box_raw.split()
        for i in range(len(virial_split)):
          virial_list.append(float(virial_split[i].decode()))
        virial = np.array(virial_list)
        virial_array = virial.reshape(frames_num, 9)
      else:
        virial_array = np.array([[]])
    else:
      print ('Need coord.raw, box.raw, force.raw, and energy.raw files, lack of essential file!')
      exit()

    return energy_array, coord_array, frc_array, box_array, virial_array

def raw_data_to_set(parts, shuffle_data, data_dir, energy_array, coord_array, frc_array, box_array, virial_array):

  '''
  raw_data_to_set: divide the raw data into several parts

  Args:
    parts: int
      parts is the number of divided parts.
    shuffle_data: bool
      shuffle_data is whether we need to shuffle data.
    data_dir: string
      data_dir is the directory of raw data.
    energy_array: 1-d float array, dim = num of frames
      energy_array is the 1-d array of energies.
    coord_array: 2-d float array, dim = (num of frames)*(3*(num of atoms))
      coord_array is the 2-d array of coordinates
    frc_array: 2-d float array, dim = (num of frames)*(3*(num of atoms))
      frc_array is the 2-d array of forces.
    box_array: 2-d float array, dim = (num of frames)*9
      box_array is the 2-d array of cells.
    virial_array: 2-d float array, dim = (num of frames)*9
      virial_array is the 2-d array of virials
  Returns :
    train_data_num : int
      train_data_num is the number of training data.
  '''

  frames_num = len(energy_array)
  index = np.array(list(range(frames_num)))
  if shuffle_data:
    np.random.shuffle(index)
  part_num = math.ceil(frames_num/parts)
  index_parts = []
  index_temp = data_op.list_split(list(index), part_num)
  for i in index_temp:
    index_parts.append(i)

  train_data_num = 0

  if ( len(index_parts) == 1 ):
    train_data_num = frames_num

  for i in range(len(index_parts)):
    if ( len(index_parts) > 1 and i != len(index_parts)-1 ):
      train_data_num = train_data_num+len(index_parts[i])
    if ( i > 10 or i == 10 ):
      sub_dir_name = 'set.0'+str(i)
    elif ( i < 10 ):
      sub_dir_name = 'set.00'+str(i)
    elif ( i > 100 ):
      sub_dir_name = 'set.'+str(i)

    cmd = "mkdir %s" % (sub_dir_name)
    call.call_simple_shell(data_dir, cmd)

    energy_array_i = energy_array[index_parts[i]].astype(np.float32)
    coord_array_i = coord_array[index_parts[i]].astype(np.float32)
    frc_array_i = frc_array[index_parts[i]].astype(np.float32)
    box_array_i = box_array[index_parts[i]].astype(np.float32)

    np.save(''.join((data_dir, '/', sub_dir_name, '/energy.npy')), energy_array_i)
    np.save(''.join((data_dir, '/', sub_dir_name, '/coord.npy')), coord_array_i)
    np.save(''.join((data_dir, '/', sub_dir_name, '/force.npy')), frc_array_i)
    np.save(''.join((data_dir, '/', sub_dir_name, '/box.npy')), box_array_i)

    #virial is not necessary
    if ( virial_array.shape != (1,0) ):
      virial_array_i = virial_array[index_parts[i]].astype(np.float32)
      np.save(''.join((data_dir, '/', sub_dir_name, '/virial.npy')), virial_array_i)

  return train_data_num

if __name__ == '__main__':
  from CP2K_kit.deepff import load_data

  #Test load_data_from_sepfile function
#  file_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir/iter_1/03.cp2k_calc/cp2k_task_1'
#  file_prefix = 'task_'
#  proj_name = 'cp2k'
#  load_data.load_data_from_sepfile(file_dir, file_prefix, proj_name)
#  load_data.raw_to_set(file_dir + '/data', 1)
  #Test load_data_from_dir and raw_to_set functions
  work_dir = '/home/lujunbo/code/github/CP2K_kit/deepff/work_dir'
  save_dir_name = 'init_train_data'
  save_dir = ''.join((work_dir, '/', save_dir_name))
  save_dir = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2_metad/iter_1/03.cp2k_calc/sys_0/data'
#  load_data.load_data_from_dir('/home/lujunbo/WORK/Deepmd/deepmd-kit/64_H2O', work_dir, save_dir, 'WATER')
  load_data.raw_to_set(save_dir, 1)
