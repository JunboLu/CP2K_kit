#! /usr/env/bin python

import os
import sys
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import data_op

traj_coord_file = str(sys.argv[1])
traj_box_file = str(sys.argv[2])
input_file = str(sys.argv[3])

atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
traj_info.get_traj_info(traj_coord_file, 'coord_xyz')

work_dir = os.getcwd()

for i in range(frames_num):
  task_dir = ''.join((work_dir, '/task_', str(i)))
  if ( not os.path.exists(task_dir) ):
    cmd = "mkdir %s" %(''.join(('task_', str(i))))
    call.call_simple_shell(work_dir, cmd)
  cmd = "cp %s %s" %(input_file, task_dir)
  call.call_simple_shell(work_dir, cmd)
  cood_file_name = ''.join((task_dir, '/coord'))
  box_file_name = ''.join((task_dir, '/box'))
  coord_file = open(cood_file_name, 'w')
  box_file = open(box_file_name, 'w')
  for j in range(atoms_num):
    line = linecache.getline(traj_coord_file, i*(atoms_num+2)+j+3)
    coord_file.write(line)
  coord_file.close()
  line = linecache.getline(traj_box_file, i+2)
  line_split = data_op.split_str(line, ' ')
  box_file.write(''.join(('A  ', line_split[2], '  ', line_split[3], '  ', line_split[4], '\n')))
  box_file.write(''.join(('B  ', line_split[5], '  ', line_split[6], '  ', line_split[7], '\n')))
  box_file.write(''.join(('C  ', line_split[8], '  ', line_split[9], '  ', line_split[10], '\n')))
  box_file.close()

linecache.clearcache()
