#!/usr/bin/env python

import os
import math
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import list_dic_op

def get_block_base(file_name):

  '''
  get_block_base : get several important information of trajectory file.

  Args :
    file_name : string
      file_name is the name of trajectory file used to analyze.
  Returns :
    b_num : int
      b_num is the number of lines in a block in trajectory file.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    base : int
      base is the number of lines before structure in a structure block.
    file_start : int
      file_start is the starting frame in trajectory file.
  '''

  if (".xyz" in file_name):
    line = linecache.getline(file_name, 1)
    b_num = int(line.strip('\n'))
    p_base = 0
    base = 2
    line = linecache.getline(file_name, 2)
    line_split = list_dic_op.str_split(line, ' ')
    if ( "-pos-" in file_name or "-vel-" in file_name or "-frc-" in file_name):
      file_start = int(line_split[2].strip(','))
    else:
      file_start = 0

  if ("-mix-1.ener" in file_name):
    b_num = 1
    p_base = 0
    base = 0
    line = linecache.getline(file_name, 1)
    line_split = list_dic_op.str_split(line, ' ')
    file_start = int(line_split[0])

  if ("-1.ener" in file_name and "mix" not in file_name):
    b_num = 1
    p_base = 1
    base = 0
    line = linecache.getline(file_name, 2)
    line_split = list_dic_op.str_split(line, ' ')
    file_start = int(line_split[0])

  if ("LagrangeMultLog" in file_name):
    b_num = 2
    p_base = 0
    base = 0
    file_start = 0

  return b_num, p_base, base, file_start

def find_breakpoint(file_name):

  '''
  find_breakpoint : find the incomplete frame in a trajectory file.

  Args :
    file_name : string
      file_name is the name of trajectory file used to analyze.
  Returns :
    breakpoint_frame : 1-d int list
  '''

  blocks_num, pre_base, base, frame_start = get_block_base(file_name)

  work_dir = os.getcwd()
  cmd = "grep %s %s" % ("'i = '", file_name)
  frames_num = len(call.call_returns_shell(work_dir, cmd))

  whole_line_num = len(open(file_name).readlines())
  breakpoint_frame = []

  #In general, we think the first frame is always correct.
  #compare_list is used for comparing. It contains atom names.
  compare_list = []
  for i in range(blocks_num):
    line = linecache.getline(file_name, pre_base+base+i+1)
    line_split = list_dic_op.str_split(line, ' ')
    compare_list.append(line_split[0])

  #This part is complicated. Maybe we need a better way.
  breakpoint = 0
  for i in range(frames_num-1):
    frame_list = []
    for j in range(blocks_num):
      line_num = pre_base+(blocks_num+base)*(i-len(breakpoint_frame))+breakpoint+base+j+1
      if ( line_num == whole_line_num ):
        breakpoint_frame.append(i+frame_start)
        break
      else:
        line = linecache.getline(file_name,line_num)
        c = list_dic_op.str_split(line, ' ')
        frame_list.append(c[0])
    if (frame_list != compare_list and frame_list != []):
      breakpoint_frame.append(i+frame_start)
      for k in range(blocks_num):
        if (frame_list[k] != compare_list[k]):
          breakpoint = breakpoint+k+base
          break

  if (breakpoint_frame == []):
    final_block = whole_line_num-(frames_num-1)*(blocks_num+base)-pre_base
    if (final_block != blocks_num+base):
      breakpoint_frame.append(frame_start+frames_num-1)

  if (breakpoint_frame != []):
    final_block = whole_line_num-(frames_num-1-len(breakpoint_frame))*(blocks_num+base)-breakpoint-pre_base
    if (final_block != blocks_num+base):
      breakpoint_frame.append(frame_start+frames_num-1)

  return breakpoint_frame

def delete_duplicate(file_name):

  '''
  delete_duplicate : delete duplicate frames in a trajectory file.

  Args :
    file_name : string
      file_name is the name of trajectory file used to analyze.
  Returns :
    none
  '''

  blocks_num, pre_base, base, frame_start = get_block_base(file_name)

  whole_line_num = len(open(file_name).readlines())
  frames_num = math.ceil((whole_line_num-pre_base)/(blocks_num+base))
  total_line = []
  #Choose the specific line
  for i in range(frames_num):
    if (".xyz" in file_name):
      line = linecache.getline(file_name, i*(blocks_num+base)+2+pre_base)
    if ("-mix-1.ener" in file_name):
      line = linecache.getline(file_name, i*(blocks_num+base)+pre_base)
    if ("-1.ener" in file_name and "mix" not in file_name):
      line = linecache.getline(file_name, i*(blocks_num+base)+pre_base)
    total_line.append(line)

  duplicate = []
  for index, line in enumerate(total_line):
    for i in range(frames_num-index-1):
      if ( line[0:13] == total_line[index+i+1][0:13] ):
        duplicate.append(index)

  #delete each duplicate frame
  #For each delete, the file will change, so (duplicate[i]-i) is the trick.
  work_dir = os.getcwd()
  for i in range(len(duplicate)):
    start_del = (duplicate[i]-i)*(blocks_num+base)+pre_base+1
    end_del = start_del+blocks_num+base-1
    cmd = "sed -ie %s %s" % (''.join((str(start_del), ',', str(end_del), 'd')), file_name)
    call.call_simple_shell(work_dir, cmd)

def choose_str(atoms_num, pre_base, base, each, start, end, file_start, file_name, choose_line, work_dir):

  '''
  choose_str : get the coordinates or velocities of choosed atoms.

  Args :
    atoms_num : int
      atoms_num is the number of atoms in trajectory file.
    pre_base : int
      pre_base is the number of lines before block of trajectory file.
    base : int
      base is the number of lines before structure in a structure block.
    each : int
      each is printing frequency of md.
    start : int
      start is the starting frame used to choose.
    end : int
      end is the ending frame used to choose.
    file_start : int
      file_start is the starting frame in trajectory file.
    file_name : string
      file_name is the name of trajectory file used to analyze.
    choose_line : 1-d int list
      choose_line is the choosed atom id.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  choose_file = ''.join((work_dir, '/choose_structure'))
  file_md = open(choose_file, 'w')

  for i in range(int((end-start)/each+1)):
    file_md.write(str(len(choose_line))+'\n')
    second_line = linecache.getline(file_name, (start-file_start+i)*(atoms_num+base)+2+pre_base)
    file_md.write(second_line)
    for j in choose_line:
      line = linecache.getline(file_name, (start-file_start+i)*(atoms_num+base)+j+pre_base+base)
      file_md.write(line)
