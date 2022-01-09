#! /bin/env/python

import os
import sys
import linecache

element = str(sys.argv[1])
file_name = str(sys.argv[2])
basis_num = len(open(file_name, 'r').readlines())
direc = os.getcwd()

line_1 = linecache.getline(file_name, 1)
line_1_split = line_1.split(' ')
a = []
for i in line_1_split:
  if ( i != '' ):
    a.append(i)

gau_num = len(a)-3

for i in range(basis_num):
  os.environ['var'] = str(i+1)
  os.system("mkdir BS_$var")

for i in range(basis_num):
  os.environ['var_1'] = str(i+1)
  os.system("cp input.inp POTENTIAL BASIS_MOLOPT BS_$var_1")

BS = []

for i in range(basis_num):
  tmp = linecache.getline(file_name, i+1)
  BS.append(tmp)

total_gau = 4*gau_num

for i in range(basis_num):
  c = []
  d = BS[i].split(' ')
  for j in range(len(d)):
    if (d[j] != ''):
      c.append(d[j])
  c[len(c)-1] = c[len(c)-1].strip('\n')
  filename = open('/%s/BS_%d/BASIS_MOLOPT' %(direc,i+1), 'a')
  filename.write('%s  UBS\n' %element)
  filename.write('%d\n' %total_gau)
  for k in range(gau_num):
    filename.write('1 0 0 1 1\n')
    filename.write('  '+c[k+3] + '  1.0\n')
  for j in range(gau_num):
    filename.write('1 1 1 1 1\n')
    filename.write('  '+c[j+3].strip(',') + '  1.0\n')
  for l in range(gau_num):
    filename.write('1 2 2 1 1\n')
    filename.write('  '+c[l+3].strip(',') + '  1.0\n')
  for m in range(gau_num):
    filename.write('1 3 3 1 1\n')
    filename.write('  '+c[m+3].strip(',') + '  1.0\n')
  filename.close()





