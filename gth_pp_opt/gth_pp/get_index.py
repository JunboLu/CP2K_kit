#! /usr/env/bin python

import linecache
from CP2K_kit.tools import data_op

def get_index(file_name):
  total_step = len(open(file_name).readlines())
  value = []
  for i in range(total_step):
    line = linecache.getline(file_name, i+1)
    b = []
    a = line.split(' ')
    for k in range(len(a)):
      if (a[k] != ''):
        b.append(a[k])
    value.append(float(b[len(b)-1]))
  min_index = value.index(min(value))
  line = linecache.getline(file_name, min_index+1)
  num = data_op.get_str_num(line[5:8])

  return num

if __name__ == '__main__':
  import sys
  file_name = str(sys.argv[1])
  m=get_index(file_name)
  print (m)
