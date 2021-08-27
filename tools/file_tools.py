#! /usr/env/bin python

import os
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
