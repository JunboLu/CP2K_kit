#! /usr/env/bin python

import linecache

for i in range(3001):
  cood_file_name = ''.join(('/home/lujunbo/WORK/deepmd/c2h6/mtd/train_4_sys/cp2k_md/sys_4/data/task_', str(i), '/coord.inc'))
  coord_file = open(cood_file_name, 'w')
  for j in range(8):
    line = linecache.getline('test-pos-1.xyz', i*(8+2)+j+3)
    coord_file.write(line)
  coord_file.close()
