#! /env/bin python

import os
import linecache
import subprocess
from random import random
from subprocess import check_call

max_cycle = 100000

def split_str(str_tmp, space_char, strip_char=''):

  str_tmp_split = str_tmp.split(space_char)
  list_tmp = []
  for i in range(len(str_tmp_split)):
    if ( str_tmp_split[i] != ''):
      list_tmp.append(str_tmp_split[i])

  if ( strip_char != '' ):
    if ( list_tmp[len(list_tmp)-1] == strip_char ):
      list_tmp.remove(list_tmp[len(list_tmp)-1])
    else:
      list_tmp[len(list_tmp)-1] = list_tmp[len(list_tmp)-1].strip(strip_char)

  return list_tmp

line_246 = linecache.getline('atom.out', 246)
line_246_split = split_str(line_246, ' ')
eigen_1_pre = float(line_246_split[5].split('[')[0])

line_249 = linecache.getline('atom.out', 249)
line_249_split = split_str(line_249, ' ')
eigen_2_pre = float(line_249_split[5].split('[')[0])

line_252 = linecache.getline('atom.out', 252)
line_252_split = split_str(line_252, ' ')
eigen_3_pre = float(line_252_split[5].split('[')[0])

line_254 = linecache.getline('atom.out', 254)
line_254_split = split_str(line_254, ' ')
eigen_4_pre = float(line_254_split[5].split('[')[0])

line_244 = linecache.getline('atom.out', 244)
line_244_split = split_str(line_244, '[')
charge_1_pre = float(split_str(line_244_split[1], ' ', '')[1])

line_250 = linecache.getline('atom.out', 250)
line_250_split = split_str(line_250, '[')
charge_2_pre = float(split_str(line_250_split[1], ' ', '')[1])

parameter_1 = []
parameter_2 = []

line_3 = linecache.getline('GTH-PARAMETER', 3)
line_3_split = split_str(line_3, ' ', '\n')
r_loc = float(line_3_split[0])
parameter_1.append(float(line_3_split[2]))
parameter_1.append(float(line_3_split[3]))

line_5 = linecache.getline('GTH-PARAMETER', 5)
line_5_split = split_str(line_5, ' ', '\n')
parameter_2.append(float(line_5_split[0]))
parameter_1.append(float(line_5_split[2]))
parameter_1.append(float(line_5_split[3]))
parameter_1.append(float(line_5_split[4]))
line_6 = linecache.getline('GTH-PARAMETER', 6)
line_6_split = split_str(line_6, ' ', '\n')
parameter_1.append(float(line_6_split[0]))
parameter_1.append(float(line_6_split[1]))
line_7 = linecache.getline('GTH-PARAMETER', 7)
line_7_split = split_str(line_7, ' ', '\n')
parameter_1.append(float(line_7_split[0]))

line_8 = linecache.getline('GTH-PARAMETER', 8)
line_8_split = split_str(line_8, ' ', '\n')
parameter_2.append(float(line_8_split[0]))
parameter_1.append(float(line_8_split[2]))
parameter_1.append(float(line_8_split[3]))
parameter_1.append(float(line_8_split[4]))
line_9 = linecache.getline('GTH-PARAMETER', 9)
line_9_split = split_str(line_9, ' ', '\n')
parameter_1.append(float(line_9_split[0]))
parameter_1.append(float(line_9_split[1]))
line_10 = linecache.getline('GTH-PARAMETER', 10)
line_10_split = split_str(line_10, ' ', '\n')
parameter_1.append(float(line_10_split[0]))

line_11 = linecache.getline('GTH-PARAMETER', 11)
line_11_split = split_str(line_11, ' ', '\n')
parameter_2.append(float(line_11_split[0]))
parameter_1.append(float(line_11_split[2]))
parameter_1.append(float(line_11_split[3]))
line_12 = linecache.getline('GTH-PARAMETER', 12)
line_12_split = split_str(line_12, ' ', '\n')
parameter_1.append(float(line_12_split[0]))

direc = os.path.abspath('.')

for i in range(max_cycle):
  randlist_1 = []
  randlist_2 = []
  for j in range(17):
    randlist_1.append(random()*0.1-0.05)
  for j in range(3):
    randlist_2.append(random()*0.002-0.001)
  C1 = parameter_1[0] + randlist_1[0]
  C2 = parameter_1[1] + randlist_1[1]
  h_s_11 = parameter_1[2] + randlist_1[2]
  h_s_12 = parameter_1[3] + randlist_1[3]
  h_s_13 = parameter_1[4] + randlist_1[4]
  h_s_22 = parameter_1[5] + randlist_1[5]
  h_s_23 = parameter_1[6] + randlist_1[6]
  h_s_33 = parameter_1[7] + randlist_1[7]
  h_p_11 = parameter_1[8] + randlist_1[8]
  h_p_12 = parameter_1[9] + randlist_1[9]
  h_p_13 = parameter_1[10] + randlist_1[10]
  h_p_22 = parameter_1[11] + randlist_1[11]
  h_p_23 = parameter_1[12] + randlist_1[12]
  h_p_33 = parameter_1[13] + randlist_1[13]
  h_d_11 = parameter_1[14] + randlist_1[14]
  h_d_12 = parameter_1[15] + randlist_1[15]
  h_d_22 = parameter_1[16] + randlist_1[16]

  r_s = parameter_2[0] + randlist_2[0]
  r_p = parameter_2[1] + randlist_2[1]
  r_d = parameter_2[2] + randlist_2[2]

  param_direc = ''.join((direc, '/', 'param_', str(i)))
  if ( os.path.exists(param_direc) ):
    pass
  else:
    cmd = "mkdir %s" % (param_direc)
    check_call(cmd, cwd=direc, shell=True)

  cmd = "cp atom.inp %s" %(param_direc)
  check_call(cmd, cwd=direc, shell=True)

  param_file_name = ''.join((param_direc, '/GTH-PARAMETER'))
  param_file = open(param_file_name, 'w')
  param_file.write('Am GTH-PBE-q11\n')
  param_file.write('    4    6    1    0\n')
  param_file.write('    %f    2    %f    %f\n' %(r_loc, C1, C2))
  param_file.write('    3\n')
  param_file.write('    %f    3    %f    %f    %f\n' %(r_s, h_s_11, h_s_12, h_s_13))
  param_file.write('                     %f    %f\n' %(h_s_22, h_s_23))
  param_file.write('                           %f\n' %(h_s_33))
  param_file.write('    %f    3    %f    %f    %f\n' %(r_p, h_p_11, h_p_12, h_p_13))
  param_file.write('                     %f    %f\n' %(h_p_22, h_p_23))
  param_file.write('                           %f\n' %(h_p_33))
  param_file.write('    %f    2    %f    %f\n' %(r_d, h_d_11, h_d_12))
  param_file.write('                     %f\n' %(h_d_22))
  param_file.close()

  inp_file = ''.join((param_direc, '/atom.inp'))
  out_file = ''.join((param_direc, '/atom.out'))
  err_file = ''.join((param_direc, '/atom.err'))
  cmd = "/home/lujunbo/bin/cp2k-6.1-sopt/cp2k-6.1-Linux-x86_64.sopt %s 1> %s 2> %s" %(inp_file, out_file, err_file)
  subprocess.run(cmd, cwd=param_direc, shell=True)

  line_246 = linecache.getline(out_file, 246)
  line_246_split = split_str(line_246, ' ')
  eigen_1 = float(line_246_split[5].split('[')[0])

  line_249 = linecache.getline(out_file, 249)
  line_249_split = split_str(line_249, ' ')
  eigen_2 = float(line_249_split[5].split('[')[0])

  line_252 = linecache.getline(out_file, 252)
  line_252_split = split_str(line_252, ' ')
  eigen_3 = float(line_252_split[5].split('[')[0])

  line_254 = linecache.getline(out_file, 254)
  line_254_split = split_str(line_254, ' ')
  eigen_4 = float(line_254_split[5].split('[')[0])

  line_244 = linecache.getline(out_file, 244)
  line_244_split = split_str(line_244, '[')
  charge_1 = float(split_str(line_244_split[1], ' ', '')[1])

  line_250 = linecache.getline(out_file, 250)
  line_250_split = split_str(line_250, '[')
  charge_2 = float(split_str(line_250_split[1], ' ', '')[1])

  if ( abs(eigen_1) < abs(eigen_1_pre) and \
       abs(eigen_2) < abs(eigen_2_pre) and \
       abs(eigen_3) < abs(eigen_3_pre) and \
       abs(eigen_4) < abs(eigen_4_pre) and \
       abs(charge_1) < abs(charge_1_pre) and \
       abs(charge_2) < abs(charge_2_pre)):
    break
