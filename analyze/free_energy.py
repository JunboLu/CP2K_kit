#!/usr/bin/env python

import csv
import linecache
from CP2K_kit.tools import list_dic_op

def ti_method(force_cmd_file):

  # The unit of force is Hartree/bohr, Hartree/rad

  total_num = len(open(force_cmd_file).readlines())

  line_1 = linecache.getline(force_cmd_file, 1)
  line_1_split = list_dic_op.str_split(line_1, ' ')
  target_1 = float(line_1_split[0])

  line_2 = linecache.getline(force_cmd_file, 2)
  line_2_split = list_dic_op.str_split(line_2, ' ')
  target_2 = float(line_2_split[0])

  increment = target_2-target_1

  target_value = []
  force_value = []
  free_energy_value = []

  for i in range(total_num):
    line = linecache.getline(force_cmd_file,i+1)
    line_split = list_dic_op.str_split(line, ' ')
    target_value.append(float(line_split[0]))
    force_value.append(float(line_split[1].strip('\n')))

  for i in range(total_num):
    sum_value = 0.0
    if (i == 0):
      sum_value = 0.0
    else:
      for j in range(i):
        sum_value = sum_value+(force_value[j]+force_value[j+1])/2*increment
    #The unit of energy is kcal/mol
    free_energy_value.append((0.0-sum_value)*627.5094*1.8897259886)

  return target_value, free_energy_value

def redox_pka_slow_growth(vertical_ene,increment):

  integral = 0.0
  for i in range(len(vertical_ene)-1):
#    integral = integral+(vertical_ene[i+1]+vertical_ene[i])/2.0*increment
    integral = integral+vertical_ene[i]*increment

  with open("mix_ene.csv","w") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["time","vertical_ene"])
    for i in range(len(vertical_ene)):
      writer.writerow([i*increment,vertical_ene[i]])

  return integral

def free_energy_run(free_energy_param, work_dir):

  if ( 'method' in free_energy_param.keys() ):
    method = free_energy_param['method']
  else:
    print ('No free energy method found, please set method')
    exit()

  if ( method == 'ti' ):
    if ( 'ti_file' in free_energy_param.keys() ):
      ti_file = free_energy_param['ti_file']
      target, free_energy_value = free_energy.cmd_method(ti_file)
      ti_free_energy_file = ''.join((work_dir, '/ti_free_energy.csv'))
      with open(ti_free_energy_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['target_value', 'free_energy'])
        for i in range(len(target)):
          writer.writerow([garget[i],free_energy_value[i]])
    else:
      print ('No TI file found, please set ti_file')
      exit()
