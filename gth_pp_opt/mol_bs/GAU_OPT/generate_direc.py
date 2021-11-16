#! /bin/env/python

import subprocess

total_exp = 25
total_range = 40
pre = [0.1,0.2,0.3,0.4,0.5,0.6]

for pre_value in pre:
  cmd = "mkdir %s" %(''.join(('pre_', str(pre_value))))
  subprocess.run(cmd, cwd='./', shell=True)
  for i in range(total_range):
    cmd = "mkdir %s" %(''.join(('pre_', str(pre_value), '/range_', str(i+1))))
    subprocess.run(cmd, cwd='./', shell=True)
    for j in range(total_exp):
      cmd = "mkdir %s" %(''.join(('pre_', str(pre_value), '/range_', str(i+1), '/exp_', str(j+1))))
      subprocess.run(cmd, cwd='./', shell=True)

for pre_value in pre:
  for i in range(total_range):
    for j in range(total_exp):
      cmd = "cp input.inp POTENTIAL %s" %(''.join(('pre_', str(pre_value), '/range_', str(i+1), '/exp_', str(j+1))))
      subprocess.run(cmd, cwd='./', shell=True)

for pre_value in pre:
  for i in range(total_range):
    for j in range(total_exp):
      file_name = ''.join(('pre_', str(pre_value), '/range_', str(i+1), '/exp_', str(j+1), '/input.inp'))
      cmd = "sed -ie '35s/.*/CONFINEMENT %f %d %d/' %s" %(pre_value, i+1, j+1, file_name)
      subprocess.run(cmd, cwd='./', shell=True)
