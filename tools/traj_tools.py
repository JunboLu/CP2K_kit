#!/usr/bin/env python

import os
import math
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op

def get_block_base(file_name, file_type):

  '''
  get_block_base: get several important information of trajectory file.

  Args:
    file_name: string
      file_name is the name of trajectory file used to analyze.
    file_type: string
      file_type is the type of file.
  Returns:
    block_num: int
      block_num is the number of lines in a block in trajectory file.
    pre_base: int
      pre_base is the number of lines before block of trajectory file.
    base: int
      base is the number of lines before structure in a structure block.
    file_start: int
      file_start is the starting frame in trajectory file.
  '''

  if ( file_type == 'coord_xyz' or file_type == 'vel' or file_type == 'frc' ):
    line = linecache.getline(file_name, 1)
    block_num = int(line.strip('\n'))
    pre_base = 0
    pre_base_block = 2
    end_base_block = 0
    line = linecache.getline(file_name, 2)
    line_split = data_op.split_str(line, ' ')
    if ( len(line_split) > 1 ):
      if ( data_op.eval_str(line_split[2].strip(',')) == 1 ):
        file_start = int(line_split[2].strip(','))
      else:
        file_start = 0
    else:
      file_start = 0

    linecache.clearcache()

  if ( file_type == 'coord_pdb' ):
    pre_base = 0
    while True:
      line = linecache.getline(file_name, pre_base+1)
      line_split = data_op.split_str(line, ' ')
      if ( line_split[0] != 'ATOM' ):
        pre_base = pre_base+1
      else:
        break
    block_num = 0
    while True:
      line = linecache.getline(file_name, pre_base+block_num+1)
      line_split = data_op.split_str(line, ' ')
      if ( line_split[0] == 'ATOM' ):
        block_num = block_num+1
      else:
        break
    end_base_block = 1
    whole_line_num = len(open(file_name, 'r').readlines())
    if ( whole_line_num > pre_base+block_num+1 ):
      pre_base_block = 0
      while True:
        line = linecache.getline(file_name, pre_base+block_num+end_base_block+pre_base_block+1)
        line_split = data_op.split_str(line, ' ')
        if ( line_split[0] != 'ATOM' ):
          pre_base_block = pre_base_block+1
        else:
          break
    else:
      pre_base_block = 0
    pre_base = pre_base-pre_base_block
    line = linecache.getline(file_name, pre_base+1)
    line_split = data_op.split_str(line, ' ')
    file_start = int(line_split[2].strip(','))

  if ( file_type == 'mix_ener' ):
    block_num = 1
    pre_base = 0
    pre_base_block = 0
    end_base_block = 0
    line = linecache.getline(file_name, 1)
    line_split = data_op.split_str(line, ' ')
    file_start = int(line_split[0])

    linecache.clearcache()

  if ( file_type == 'ener' ):
    block_num = 1
    pre_base = 1
    pre_base_block = 0
    end_base_block = 0
    line = linecache.getline(file_name, 2)
    line_split = data_op.split_str(line, ' ')
    file_start = int(line_split[0])

    linecache.clearcache()

  if ( file_type == 'lagrange' ):
    block_num = 2
    pre_base = 0
    pre_base_block = 0
    end_base_block = 0
    file_start = 0

  return block_num, pre_base, pre_base_block, end_base_block, file_start

def find_breakpoint(file_name, file_type):

  '''
  find_breakpoint: find the incomplete frame in a trajectory file.

  Args:
    file_name: string
      file_name is the name of trajectory file used to analyze.
    file_type: string
      file_type is the type of file.
  Returns:
    breakpoint: int
      breakpoint is the incomplete frame id.
  '''

  blocks_num, pre_base, pre_base_block, end_base_block, frame_start = get_block_base(file_name, file_type)

  work_dir = os.getcwd()
  cmd = "grep %s %s" % ("'i = '", file_name)
  frames_num = len(call.call_returns_shell(work_dir, cmd))

  whole_line_num = len(open(file_name).readlines())

  #In general, we think the first frame is always correct.
  #compare_list is used for comparing. It contains atom names.
  compare_list = []
  for i in range(blocks_num):
    line = linecache.getline(file_name, pre_base+pre_base_block+i+1)
    line_split = data_op.split_str(line, ' ')
    compare_list.append(line_split[0])

  breakpoint = frame_start
  for i in range(frames_num-1):
    breakpoint = i+frame_start
    frame_list = []
    for j in range(blocks_num):
      line_num = pre_base+(pre_base_block+blocks_num+end_base_block)*i+pre_base_block+j+1
      line = linecache.getline(file_name,line_num)
      c = data_op.split_str(line, ' ')
      frame_list.append(c[0])
    if (frame_list != compare_list and frame_list != []):
      break

  linecache.clearcache()

  return breakpoint

def delete_duplicate(file_name, file_type):

  '''
  delete_duplicate: delete duplicate frames in a trajectory file.

  Args:
    file_name: string
      file_name is the name of trajectory file used to analyze.
    file_type: string
      file_type is the type of file.
  Returns:
    none
  '''

  blocks_num, pre_base, pre_base_block, end_base_block, frame_start = get_block_base(file_name, file_type)

  whole_line_num = len(open(file_name).readlines())
  frames_num = math.ceil((whole_line_num-pre_base)/(pre_base_block+blocks_num+end_base_block))
  total_line = []
  #Choose the specific line
  for i in range(frames_num):
    if ( file_type == 'coord_xyz' or file_type == 'vel' or file_type == 'frc' ):
      line = linecache.getline(file_name, i*(pre_base_block+blocks_num+end_base_block)+2+pre_base)
    if ( file_type == 'mix_ener' or file_type == 'ener' ):
      line = linecache.getline(file_name, i*(pre_base_block+blocks_num+end_base_block)+pre_base+1)
    total_line.append(line)

  linecache.clearcache()

  duplicate = []
  for index, line in enumerate(total_line):
    for i in range(frames_num-index-1):
      if ( line[0:13] == total_line[index+i+1][0:13] ):
        duplicate.append(index)

  #delete each duplicate frame
  #For each delete, the file will change, so (duplicate[i]-i) is the trick.
  work_dir = os.getcwd()
  for i in range(len(duplicate)):
    start_del = (duplicate[i]-i)*(pre_base_block+blocks_num+end_base_block)+pre_base+1
    end_del = start_del+pre_base_block+blocks_num+end_base_block-1
    cmd = "sed -ie %s %s" % (''.join((str(start_del), ',', str(end_del), 'd')), file_name)
    call.call_simple_shell(work_dir, cmd)

def choose_str(atoms_num, pre_base, pre_base_block, end_base_block, each, init_step, end_step, \
               start_frame_id, traj_coord_file, choose_line, work_dir, choose_file_name):

  '''
  choose_str: get the coordinates or velocities of choosed atoms.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    pre_base: int
      pre_base is the number of lines before block of the trajectory.
    base: int
      base is the number of lines before structure in a structure block.
    each: int
      each is printing frequency of md.
    init_step: int
      init_step is the initial step frame id.
    end_step: int
      end_step is the endding step frame id.
    start_frame_id: int
      start_frame_id is the starting frame id in trajectory file.
    traj_coord_file: string
      traj_coord_file is the name of coordination trajectory file.
    choose_line: 2-d int list
      choose_line is the choosed atom id.
    work_dir: string
      work_dir is working directory of CP2K_kit.
    choose_file_name: string
      choose_file_name is the name of generated file.
  Returns:
    choose_file: string
      choose_file is the name of choosed structure.
  '''

  choose_file = ''.join((work_dir, '/', choose_file_name))
  file_md = open(choose_file, 'w')

  choose_num = []
  for i in range(len(choose_line)):
    choose_num.append(len(choose_line[i]))

  for i in range(int((end_step-init_step)/each+1)):
    file_md.write(str(sum(choose_num))+'\n')
    second_line = linecache.getline(traj_coord_file, (init_step-start_frame_id+i)*(pre_base_block+atoms_num+end_base_block)+2+pre_base)
    file_md.write(second_line)
    for j in range(len(choose_line)):
      for k in choose_line[j]:
        line = linecache.getline(traj_coord_file, (init_step-start_frame_id+i)*(pre_base_block+atoms_num+end_base_block)+k+pre_base+pre_base_block)
        file_md.write(line)

  linecache.clearcache()

  return choose_file

def order_traj_file(atoms_num, frames_num, each, init_step, traj_file, file_type, order_list, work_dir, file_name):

  '''
  order_traj_file : reorder the trajectory file.

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    frames_num: int
      frames_num is the number of frames in the trajectory file.
    each: int
      each is printing frequency of md.
    init_step: int
      init_step is the initial step frame id.
    traj_file: string
      traj_file is the name of trajectory file.
    file_type: string
      file_type is the type of file.
    order_list: 2-d int list
      order_list is the order list.
    work_dir: string
      work_dir is working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    new_traj_file_name: string
      new_traj_file_name is the name of ordered structure.
  '''

  block_num, pre_base, pre_base_block, end_base_block, start_frame_id = get_block_base(traj_file, file_type)
  new_traj_file_name = ''.join((work_dir, '/', file_name))
  new_traj_file = open(new_traj_file_name, 'w')

  for i in range(frames_num):
    line_i_1 = linecache.getline(traj_file, (pre_base_block+atoms_num+end_base_block)*int((init_step-start_frame_id)/each+i)+1+pre_base)
    line_i_2 = linecache.getline(traj_file, (pre_base_block+atoms_num+end_base_block)*int((init_step-start_frame_id)/each+i)+2+pre_base)
    new_traj_file.write(line_i_1)
    new_traj_file.write(line_i_2)
    for j in range(len(order_list)):
      for k in order_list[j]:
        line_ijk = linecache.getline(traj_file, (pre_base_block+atoms_num+end_base_block)*int((init_step-start_frame_id)/each+i)+pre_base_block+k+pre_base)
        new_traj_file.write(line_ijk)

  linecache.clearcache()
  new_traj_file.close()

  return new_traj_file_name

if __name__ == '__main__':
  from CP2K_kit.tools import traj_tools

  a, b, c, d, e = traj_tools.get_block_base('./pdb', 'coord_pdb')
  print (a, b, c, d, e)
#  file_name = '/home/lujunbo/WORK/test/TEST_CP2K_KIT/TEST_ANALYZE/center/traj_order.xyz'
#  breakpoint = traj_tools.find_breakpoint(file_name, 'coord_xyz')
#  print (breakpoint)
