#! /usr/env/bin python

import os
import sys
import linecache
import subprocess
from CP2K_kit.tools import call
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import file_tools
from CP2K_kit.tools import revise_cp2k_inp

##############################################
traj_coord_file = str(sys.argv[1])
traj_box_file = str(sys.argv[2])
each_user = int(sys.argv[3])
input_file = str(sys.argv[4])
##############################################

work_dir = os.getcwd()
subprocess.run('mkdir data', cwd=work_dir, shell=True)
data_dir = ''.join((work_dir, '/data'))

revise_cp2k_inp.delete_line("'WFN_RESTART_FILE_NAME'", input_file, work_dir)
cp2k_inp_file_upper = file_tools.upper_file(input_file, work_dir)
pot_line_num = file_tools.grep_line_num("'POTENTIAL_FILE_NAME'", cp2k_inp_file_upper, work_dir)[0]
call.call_simple_shell(work_dir, 'rm %s' %(cp2k_inp_file_upper))
proj_name = revise_cp2k_inp.get_proj_name(input_file, work_dir)

atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each_traj, start_frame_id, end_frame_id, time_step = \
traj_info.get_traj_info(traj_coord_file, 'coord_xyz')

each = each_user*each_traj
frames_num_stat = int((end_frame_id-start_frame_id)/each+1)

for i in range(frames_num):
  task_dir = ''.join((data_dir, '/task_', str(i)))
  if ( not os.path.exists(task_dir) ):
    cmd = "mkdir %s" %(''.join(('task_', str(i))))
    call.call_simple_shell(data_dir, cmd)
  cmd = "cp %s %s" %(input_file, task_dir)
  call.call_simple_shell(work_dir, cmd)
  inp_file = ''.join((task_dir, '/', input_file))
  if ( i != 0 ):
    cmd = "sed -i '%d s/^/    WFN_RESTART_FILE_NAME ..\/traj_%d\/%s-RESTART.wfn\\n/' %s" \
          %(pot_line_num+1, i-1, proj_name, inp_file)
    call.call_simple_shell(work_dir, cmd)

  cood_file_name = ''.join((task_dir, '/coord'))
  box_file_name = ''.join((task_dir, '/box'))
  coord_file = open(cood_file_name, 'w')
  box_file = open(box_file_name, 'w')
  for j in range(atoms_num):
    line = linecache.getline(traj_coord_file, i*(atoms_num+2)*each+j+3)
    coord_file.write(line)
  coord_file.close()
  line = linecache.getline(traj_box_file, i+2)
  line_split = data_op.split_str(line, ' ')
  box_file.write(''.join(('A  ', line_split[2], '  ', line_split[3], '  ', line_split[4], '\n')))
  box_file.write(''.join(('B  ', line_split[5], '  ', line_split[6], '  ', line_split[7], '\n')))
  box_file.write(''.join(('C  ', line_split[8], '  ', line_split[9], '  ', line_split[10], '\n')))
  box_file.close()

linecache.clearcache()
