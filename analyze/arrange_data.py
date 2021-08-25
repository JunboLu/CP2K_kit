#!/usr/bin/env python

import os
import csv
import math
import linecache
import numpy as np
from CP2K_kit.tools import call
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import traj_tools
from CP2K_kit.tools import numeric
from CP2K_kit.lib import statistic_mod
from CP2K_kit.analyze import free_energy
from CP2K_kit.analyze import check_analyze

mulliken_pre_base = 5
mulliken_late_base = 3

def arrange_temp(frames_num, pre_base, time_step, file_name, work_dir, each=1):

  '''
  arrange_temp : arrange temperatures from trajectory file.

  Args :
    frames_num : int
      frames_num is the number of frames in trajectory file.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    each : int
      each is printing frequency of md.
  Returns :
    none
  '''

  time = []
  temp = []

  for i in range(frames_num):
    time.append(time_step*i*each)
    line_i = linecache.getline(file_name, i+pre_base+1)
    line_i_split = data_op.str_split(line_i, ' ')
    temp.append(float(line_i_split[3])) #The temperature is in 4th row in energy file.

  temp_file = ''.join((work_dir, '/temperature.csv'))
  with open(temp_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'temperature'])
    for i in range(frames_num):
      writer.writerow([time[i], temp[i]])

  return temp_file

def arrange_pot(frames_num, pre_base, time_step, file_name, work_dir, each=1):

  '''
  arrange_pot : Arrange potential energies from trajectory file.

  Args :
    frames_num : int
      frames_num is the number of frames in trajectory file.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    each : int
      each is printing frequency of md.
  Returns :
    none
  '''

  time = []
  pot = []

  for i in range(frames_num):
    time.append(time_step*i*each)
    line_i = linecache.getline(file_name, i+pre_base+1)
    line_i_split = data_op.str_split(line_i, ' ')
    pot.append(float(line_i_split[4])) #The potential energy is in 5th row in energy file. 

  pot_file = ''.join((work_dir, '/potential.csv')) #The energy unit is Hartree.
  with open(pot_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'energy'])
    for i in range(frames_num):
      writer.writerow([time[i], pot[i]])

  return pot_file

def arrange_mulliken(frames_num, atoms_num, time_step, atom_id, file_name, work_dir, each=1):

  '''
  arrange_mulliken : arrange mulliken charge from trajectory file.

  Args :
    frames_num : int
      frames_num is the number of frames in trajectory file.
    atoms_num : int
      atoms_num is the number of atoms in trajectory file.
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    atom_id : int list
      atom_id is the id of atoms to be analyzed.
      Example : [1,2,3,7,8]
    file_name : string
      file_name is the name of file used to analyze.
    work_dir : string
      work_dir is working directory of CP2K_kit.
    each : int
      each is printing frequency of md.
  Returns :
    none
  '''

  time = []
  mulliken = []

  for i in range(frames_num):
    time.append(i*time_step*each)
    if ( isinstance(atom_id, int) ):
      line_i = linecache.getline(file_name, i*(mulliken_pre_base+atoms_num+mulliken_late_base)+atom_id+mulliken_pre_base)
      line_i_split = data_op.str_split(line_i, ' ')
      mulliken.append(float(line_i_split[6].strip('\n'))) #Mulliken charge is in 6th row of mulliken file.
    elif ( isinstance(atom_id, list) ):
      mulliken_i = 0.0
      for j in range(len(atom_id)):
        line_ij = linecache.getline(file_name, i*(mulliken_pre_base+atoms_num+mulliken_late_base)+atom_id[j]+mulliken_pre_base)
        line_ij_split = data_op.str_split(line_ij, ' ')
        mulliken_i = mulliken_i + float(line_ij_split[6].strip('\n'))
      mulliken.append(mulliken_i)

  mulliken_file = ''.join((work_dir, '/mulliken.csv'))
  with open(mulliken_file, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['time', 'mulliken'])
    for i in range(frames_num):
      writer.writerow([time[i], mulliken[i]])

  return mulliken_file

def arrange_vertical_energy(time_step, final_time_unit, start, end, file_start, mix_ene_file, \
                            row_ox, row_red, redox_type, slow_growth, work_dir, each=1):

  '''
  arrange_vertical_energy : arrange vertical energy from trajectory file.

  Args :
    time_step : float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    final_time_unit : string
      final_time_unit is the unit of time in the result.
    start : int
      start is the starting frame used to analyze.
    end : int
      end is the ending frame used to analyze.
    file_start : int
      file_start is the starting frame in trajectory file.
    mix_ene_file : string
      mix_ene_file is the mixing energy trajctory file.
    row_ox : int
      row_ox is the row of oxidation species in trajectory file.
    row_red : int
      row_red is the row of reduced species in trajectory file.
    redox_type : string
      redox_type is the type of redox. Two choices: oxidation, reduction.
    slow_growth : int
      slow_growth is the keyword show whether user use slow-growth method.
      Example : 0 means no slow-growth method, 1 means slow-growth method
    work_dir : string
      work_dir is working directory of CP2K_kit.
    each : int
      each is printing frequency of md.
  Returns :
    delta_e_avg : float
      delta_e_avg is the averaged value of vertical energy.
    rmse : float
      rmse is the statistical error of delta_e_avg
    vertical_ene_array : 1d float array
      vertical_ene_array is vertical energy along MD.
  '''

  vertical_ene = []
  mix_ene = []
  time = []

  #index_1 and index_2 are different when we treat different redox_type.
  if (redox_type == 'oxidation'):
    index_1 = row_ox-1
    index_2 = row_red-1
  elif (redox_type == 'reduction'):
    index_1 = row_red-1
    index_2 = row_ox-1

  stat_num = end-start+1
  if ( final_time_unit == 'fs' ):
    for i in range(stat_num):
      time.append(time_step*i*each)
  elif ( final_time_unit == 'ps' ):
    for i in range(stat_num):
      time.append(time_step*i*each*0.001)
  elif ( final_time_unit == 'ns' ):
    for i in range(stat_num):
      time.append(time_step*i*each*0.000001)

  for i in range(stat_num):
    line_i = linecache.getline(mix_ene_file, start-file_start+i+1)
    line_i_split = data_op.str_split(line_i, ' ')
    delta_ene = (float(line_i_split[index_1])-float(line_i_split[index_2]))*27.2114 #The unit for vertical energy is eV.
    vertical_ene.append(delta_ene)
    if ( slow_growth == 0 ):
      mix_ene.append(float(line_i_split[2]))

  vertical_ene_array = np.asfortranarray(vertical_ene, dtype='float32')

  #If we do not use slow-growth, we arrange data for one lamda.
  #If we use slow-growth, we arrange data and then return to calculate free energy directly.
  if ( slow_growth == 0 ):
    mix_ene_file = ''.join((work_dir, '/vertical_ene.csv'))
    with open(mix_ene_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time', 'vertical_ene', 'mix_ene'])
      for i in range(stat_num):
        writer.writerow([time[i], vertical_ene[i], mix_ene[i]])

    #Get averaged values of square of vertical energy and vertical energy.
    vertical_ene_sqr = np.array(vertical_ene)**2  #vertical_ene_sqr is the squre of vertical_ene
    vertical_ene_sqr_array = np.asfortranarray(vertical_ene_sqr, dtype='float32')
    delta_e_sqr_avg, sigma_1 = statistic_mod.statistic.numerical_average(vertical_ene_sqr_array, stat_num)
    delta_e_avg, sigma_2 = statistic_mod.statistic.numerical_average(vertical_ene_array, stat_num)
    rmse = np.sqrt(delta_e_sqr_avg - delta_e_avg**2)

    #Read csv and get vertical energy for every time
    vertical_ene_avg = []
    for i in range(len(vertical_ene)):
      sum_value = 0.0
      for j in range(i+1):
        sum_value = sum_value + vertical_ene[j]
      vertical_ene_avg.append(sum_value/(i+1))

    vertical_ene_avg_file = ''.join((work_dir, '/vertical_ene_avg.csv'))
    with open(vertical_ene_avg_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['time', 'vertical_ene_avg'])
      for i in range(len(vertical_ene)):
        writer.writerow([time[i], vertical_ene_avg[i]])

    max_vertical_ene = max(vertical_ene)
    min_vertical_ene = min(vertical_ene)
    increment = (max_vertical_ene-min_vertical_ene)/500.0
    vertical_freq = []
    for i in range(500):
      number = 0
      for j in range(len(vertical_ene)):
        if ( vertical_ene[j] > (min_vertical_ene+i*increment) and vertical_ene[j] < (min_vertical_ene+(i+1)*increment) ):
          number = number + 1
      vertical_freq.append(number/len(vertical_ene))
    vertical_freq_fit = numeric.savitzky_golay(np.array(vertical_freq), 201, 3)

    frequency_file = ''.join((work_dir, '/frequency.csv'))
    with open(frequency_file, 'w') as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(['vertical_energy','frequency','frequency_fit'])
      for i in range(len(vertical_freq)):
        writer.writerow([min_vertical_ene+i*increment, vertical_freq[i], vertical_freq_fit[i]])

    return delta_e_avg, rmse
  elif ( slow_growth == 1 ):
    return vertical_ene_array

def arrange_ti_force(stat_num, lagrange_file):

  '''
  arrange_ti_force : arrange force in thermodynamic integration calculation.

  Args :
    stat_num : int
      stat_num is the number of data used to be analyzed.
    lagrange_file : string
      lagrange_file is the file containg lagrange force.
  Returns :
    s_avg : float
      s_avg is the averaged value of lagrange force.
    sq_avg : float
      sq_avg is the averaged value of square of lagrange force.
  '''

  #Get trajectory information of lagrange file.
  blocks_num, pre_base, base, frame_start = traj_tools.get_block_base(lagrange_file, 'lagrange')
  whole_line_num = len(open(lagrange_file).readlines())
  frames_num = int((whole_line_num-pre_base)/(blocks_num+base))

  s_sum = 0.0
  sq_sum = 0.0

  for i in range(stat_num):
    line_i = linecache.getline(lagrange_file, (frames_num-stat_num+i)*(blocks_num+base)+1)
    line_i_split = data_op.str_split(line_i, ' ')
    s_sum = s_sum+float(line_i_split[3].strip('\n'))

    sq_sum = sq_sum+(float(line_i_split[3].strip('\n')))**2

  s_avg = s_sum/stat_num
  sq_avg = np.sqrt(sq_sum/stat_num-s_sum**2/stat_num**2)/(np.sqrt(stat_num))

  return s_avg, sq_avg

def arrange_data_run(arrange_data_param, work_dir):

  '''
  arrange_data_run : kernel function to run arrange_data.

  Args :
    arrange_data_param : dictionary
      arrange_data_param contains information of arrange_data.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  '''

  arrange_data_param = check_analyze.check_arrange_data_inp(arrange_data_param)

  print ('ARRANGE_DATA'.center(80, '*'), flush=True)
  #arrange temperature
  if ( 'temperature' in arrange_data_param ):
    temp_param = arrange_data_param['temperature']
    traj_ener_file = temp_param['traj_ener_file']
    blocks_num, base, pre_base, frames_num, each, start_id, end_id, time_step = \
    traj_info.get_traj_info(traj_ener_file, 'ener')

    str_print = 'Analyze temperature vs time from %s' %(traj_ener_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    temp_file = arrange_temp(frames_num, pre_base, time_step, traj_ener_file, work_dir, each)

    str_print = 'Temperature vs time is written in %s' %(temp_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

  #arrange potential energy
  elif ( 'potential' in arrange_data_param ):
    pot_param = arrange_data_param['potential']
    traj_ener_file = pot_param['traj_ener_file']
    blocks_num, base, pre_base, frames_num, each, start_id, end_id, time_step = \
    traj_info.get_traj_info(traj_ener_file, 'ener')

    str_print = 'Analyze potential energy vs time from %s' %(traj_ener_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    pot_file = arrange_pot(frames_num, pre_base, time_step, traj_ener_file, work_dir, each)

    str_print = 'Potential energy vs time is written in %s' %(pot_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

  #arrange mulliken charge
  elif ( 'mulliken' in arrange_data_param ):
    mulliken_param = arrange_data_param['mulliken']
    traj_mul_file = mulliken_param['traj_mul_file']
    atom_id = mulliken_param['atom_id']
    time_step = mulliken_param['time_step']
    each = mulliken_param['each']

    str_print = 'Analyze mulliken charge vs time from %s' %(traj_mul_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    #line 322-330 used to get the number of atoms.
    cmd = "grep -n '#  Atom  Element' %s" % (traj_mul_file)
    cmd_return = call.call_returns_shell(work_dir, cmd)
    line_num_1 = int(cmd_return[0].split(':')[0])

    cmd = "grep -n '# Total charge' %s" % (traj_mul_file)
    cmd_return = call.call_returns_shell(work_dir, cmd)
    line_num_2 = int(cmd_return[0].split(':')[0])

    atoms_num = line_num_2-line_num_1-1
    #'#  Atom  Element' and '# Total charge' are two keywords.

    whole_line_num_1 = len(open(traj_mul_file).readlines())
    frames_num = math.ceil(whole_line_num_1/(mulliken_pre_base+atoms_num+mulliken_late_base))

    mulliken_file = arrange_mulliken(frames_num, atoms_num, time_step, atom_id, traj_mul_file, work_dir, each)

    str_print = 'Mulliken charge vs time is written in %s' %(mulliken_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

  #arrange vertical energy
  elif ( 'vertical_energy' in arrange_data_param ):
    vert_ene_param = arrange_data_param['vertical_energy']

    traj_mix_ener_file = vert_ene_param['traj_mix_ener_file']
    blocks_num, base, pre_base, frames_num, each, start_id, end_id, time_step = \
    traj_info.get_traj_info(traj_mix_ene_file, 'ener')
    row_ox = vert_ene_param['row_ox']
    row_red = vert_ene_param['row_red']
    redox_type = vert_ene_param['redox_type']
    slow_growth = vert_ene_param['slow_growth']
    init_step = vert_ene_param['init_step']
    end_step = vert_ene_param['end_step']
    final_time_unit = vert_ene_param['final_time_unit']

    if ( slow_growth == 0 ):
      delta_ene, rmse = arrange_vertical_energy(time_step, final_time_unit, init_step, end_step,start_id, \
                                                traj_mix_ener_file, row_ox, row_red, redox_type, slow_growth, work_dir)
      print ('Average vertical energy is %f eV, and error is %f eV' % (delta_ene, rmse), flush=True)
    elif ( slow_growth == 1 ):
      vert_ene = arrange_vertical_energy(time_step, final_time_unit, init_step, end_step, start_id, \
                                         traj_mix_ener_file, row_ox, row_red, redox_type, slow_growth, work_dir)
      increment = vert_ene_param['increment']
      redox_pka_free_ene = free_energy.redox_pka_slow_growth(vert_ene, increment)
      print ('The redox free energy is %f ev' %(redox_pka_free_ene), flush=True)

  #arrange thermodynamic integration force
  elif ( 'ti_force' in arrange_data_param ):
    ti_force_param = arrange_data_param['ti_force']

    traj_lag_file = ti_force_param['traj_lag_file']
    stat_num = ti_force_param['stat_num']

    str_print = 'Extract lagrange force from %s' %(traj_lag_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    force_avg, error_avg = arrange_ti_force(stat_num, traj_lag_file)
    print ("The averaged force is %f and averaged error is %f" %(force_avg, error_avg), flush=True)
