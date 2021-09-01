#! /usr/env/bin python

import os
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op

def upper_file(file_name, save_dir):

  '''
  upper_file: make the character in the file to be upper

  Args:
    file_name: string
      file_name is the name of the file needed to be revised.
    save_dir: string
      save_dir is the saving dir for revised file.
  Returns:
    rev_file_name: string
      rev_file_name is the name of revised file.
  '''

  rev_file_name = ''.join((os.path.abspath(save_dir), '/upper_file'))
  rev_file = open(rev_file_name, 'w')
  origin_file = open(file_name, 'r')

  while True:
    line = origin_file.readline()
    rev_file.write(line.upper())
    if not line:
      break

  rev_file.close()

  return rev_file_name

def space_file(file_name, space_char, save_dir):

  '''
  space_file

  Args:
    file_name: string
      file_name is the name of the file needed to be revised.
    space_char: string
    save_dir: string
      save_dir is the saving dir for revised file.
  Returns:
    rev_file_name: string
      rev_file_name is the name of revised file.
  '''

  rev_file_name = ''.join((os.path.abspath(save_dir), '/space_file'))
  rev_file = open(rev_file_name, 'w')
  origin_file = open(file_name, 'r')

  while True:
    line = origin_file.readline()
    line_split = data_op.str_split(line, ' ')
    line_comb = data_op.comb_list_2_str(line_split, ' ')
    rev_file.write(line_comb)
    if not line:
      break

  rev_file.close()

  return rev_file_name

def grep_line_num(choosed_str, file_name, work_dir):

  '''
  grep_line_num : get the line number for a choosed string.

  Args:
    choosed_str: string
      choosed_str is the choosed string.
    file_name: string
      file_name is the name of file.
  Returns:
    line_num : int
      line_num is the line number for a choosed string.
  '''

  line_num = []
  cmd = "grep -n %s %s" %(choosed_str, file_name)
  line = call.call_returns_shell(work_dir, cmd)
  if ( len(line) != 0 ):
    for i in range(len(line)):
      line_num.append(int(line[i].split(':')[0]))
  else:
    line_num = 0

  return line_num

if __name__ == '__main__':

  choosed_line = "'&force_eval'"
  file_name = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md/input.inp'
  work_dir = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md'
  line_num = grep_line_num(choosed_line, file_name, work_dir)
  print (line_num)

