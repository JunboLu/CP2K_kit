#! /usr/bin/env python

import sys
import linecache
import numpy as np

file_name = str(sys.argv[1])
whole_line_num = len(open(file_name, 'r').readlines())

small_gau = []
for i in range(whole_line_num):
  c = []
  line = linecache.getline(file_name, i+1)
  d = line.split(' ')
  for j in range(len(d)):
    if (d[j] != ''):
      c.append(d[j])
  small_gau.append(float(c[3]))

new_file = open('kk_sort', 'w')
list_order_index = np.argsort(-np.array(small_gau))
for i in list_order_index:
  line = linecache.getline(file_name, i+1)
  new_file.write(line)

new_file.close()
