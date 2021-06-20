#! /bin/env/python

import os

total_exp = 25
total_range = 40
pre = [0.1,0.2,0.3,0.4,0.5,0.6]

for pre_value in pre:
  os.environ['pre_direc'] = str(pre_value)
  os.system("mkdir pre_$pre_direc")
  for i in range(total_range):
    os.environ['var'] = str(i+1)
    os.system("mkdir pre_$pre_direc/range_$var")

  for i in range(total_range):
    os.environ['var'] = str(i+1)
    for j in range(total_exp):
      os.environ['var_1'] = str(j+1)
      os.system("mkdir pre_$pre_direc/range_$var/exp_$var_1")

for pre_value in pre:
  os.environ['pre_direc'] = str(pre_value)
  for i in range(total_range):
    os.environ['var'] = str(i+1)
    for j in range(total_exp):
      os.environ['var_1'] = str(j+1)
      os.system("cp input.inp POTENTIAL pre_$pre_direc/range_$var/exp_$var_1")

  for i in range(total_range):
    os.environ['var'] = str(i+1)
    for j in range(total_exp):
      os.environ['var_1'] = str(j+1)
      os.environ['var_2'] = 'CONFINEMENT'
      os.environ['var_3'] = str(pre)
      os.environ['var_4'] = str(i+1)
      os.environ['var_5'] = str(j+1)
      os.system("sed -ie '35s/.*/'$var_2' '$var_3' '$var_4' '$var_5'/' pre_$pre_direc/range_$var/exp_$var_1/input.inp")
      os.system("rm pre_$pre_direc/range_$var/exp_$var_1/input.inpe")
