#! /bin/env/python

import os
import sys
import linecache

#basis_num = int(sys.argv[1])
#gau_num = int(sys.argv[2])
#element = sys.argv[3]
#direc = sys.argv[4]

basis_num = 10
gau_num = 7
element = 'La'
direc = os.getcwd()

for i in range(basis_num):
  os.environ['var'] = str(i+1)
  os.system("mkdir BS_$var")

for i in range(basis_num):
  os.environ['var_1'] = str(i+1)
  os.system("cp input.inp POTENTIAL cp2k.sub BASIS_MOLOPT BS_$var_1")

BS = []

for i in range(basis_num):
  tmp = linecache.getline('BASIS', i+1)
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
  filename.write('%s  CBS\n' %element)
  filename.write('%d\n' %total_gau)
  for k in range(gau_num):
    filename.write('1 0 0 1 1\n')
    filename.write('  '+c[k+4] + '  1.0\n')
  for j in range(gau_num):
    filename.write('1 1 1 1 1\n')
    filename.write('  '+c[j+4] + '  1.0\n')
  for l in range(gau_num):
    filename.write('1 2 2 1 1\n')
    filename.write('  '+c[l+4] + '  1.0\n')
  for m in range(gau_num):
    filename.write('1 3 3 1 1\n')
    filename.write('  '+c[m+4] + '  1.0\n')

  filename.close()





