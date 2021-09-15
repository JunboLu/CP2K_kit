#!/usr/bin/env python

import os
import sys
import math
import linecache
from CP2K_kit.tools import atom
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_tools

def get_traj_info(file_name, file_type, group=[[]], atom_id=[[]], return_group=False):

  '''
  get_traj_info: get several important information of trajectory

  Args:
    file_name: string
      file_name is the name of trajectory file used to analyze.
    file_type: string
      file_type is the type of file.
    group: 2-d string list
      group contain the basic atom info in a set of atoms.
      Example: [['O', 'H', 'H']]
    atom_id: 2-d int list
      atom_id is the atom id of the group atoms.
      Example: [[1,2,3,4,5,6...192]]
    return_group: bool
  Returns:
    blocks_num: int
      blocks_num is the number of lines in one block in trajectory file.
    pre_base_block: int
      pre_base_block is the number of lines before structure in a structure block.
    pre_base: int
      pre_base is the number of lines before block of trajectory file.
    frames_num: int
      frames_num is the number of frames in trajectory file.
    each: int
      each is printing frequency of md.
    start_frame_id: int
      start_frame_id is the starting frame used to choose.
    end_frame_id: int
      end_frame_id is the ending frame used to choose.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
    group_atom_1_id: 2-d int list
      group_atom_1_id is the id of first atoms in the molecules in the group.
    group_atoms_mass: 2-d float list
      group_atoms_mass contains the atoms mass for each group.
  '''

  blocks_num, pre_base, pre_base_block, end_base_block, frame_start = traj_tools.get_block_base(file_name, file_type)

  whole_line_num_1 = len(open(file_name).readlines())

  if ((whole_line_num_1-pre_base)%(pre_base_block+blocks_num+end_base_block) != 0):
    break_frame = traj_tools.find_breakpoint(file_name, file_type)
    log_info.log_error('There is incomplete frame in %s. The incomplete frame id is %d.' %(file_name, break_frame))
    exit()
  else:
    frames_num_1 = int((whole_line_num_1-pre_base)/(pre_base_block+blocks_num+end_base_block))

  if ( file_type == 'coord_xyz' or file_type == 'coord_pdb' or file_type == 'vel' or file_type == 'frc' ):
    if ( file_type == 'coord_pdb' ):
      a = linecache.getline(file_name, pre_base+1)
    else:
      a = linecache.getline(file_name, pre_base+2)
    b = data_op.split_str(a, ' ')
    if ( len(b) > 5 ):
      start_frame_id = int(b[2].strip(','))
      start_time = float(b[5].strip(','))
    else:
      start_frame_id = 0
      start_frame_time = 0.0

    if ( whole_line_num_1 > pre_base_block+blocks_num+end_base_block+pre_base ):
      if ( file_type == 'coord_pdb' ):
        a = linecache.getline(file_name, (pre_base_block+blocks_num+end_base_block)*1+pre_base+1)
      else:
        a = linecache.getline(file_name, (pre_base_block+blocks_num+end_base_block)*1+pre_base+2)
      b = data_op.split_str(a, ' ')
      if ( len(b) > 5 ):
        second_frame_id = int(b[2].strip(','))
        second_time = float(b[5].strip(','))
      else:
        second_frame_id = 0
        second_time = 0.0

      if ( file_type == 'coord_pdb' ):
        a = linecache.getline(file_name, (frames_num_1-1)*(pre_base_block+blocks_num+end_base_block)+pre_base+1)
      else:
        a = linecache.getline(file_name, (frames_num_1-1)*(pre_base_block+blocks_num+end_base_block)+pre_base+2)
      b = data_op.split_str(a, ' ')
      if ( len(b) > 5 ):
        end_frame_id = int(b[2].strip(','))
    else:
      start_frame_id = 0
      end_frame_id = 0

  if ( file_type == 'ener' or file_type == 'mix_ener' ):

    a = linecache.getline(file_name, pre_base+1)
    b = data_op.split_str(a, ' ')
    start_frame_id = int(b[0])
    start_time = float(b[1])

    if ( whole_line_num_1 > pre_base_block+blocks_num+end_base_block+pre_base+1 ):
      a = linecache.getline(file_name, (pre_base_block+blocks_num+end_base_block)*1+pre_base+1)
      b = data_op.split_str(a, ' ')
      second_frame_id = int(b[0])
      second_time = float(b[1])

      a = linecache.getline(file_name, whole_line_num_1)
      b = data_op.split_str(a, ' ')
      end_frame_id = int(b[0])
    else:
      end_frame_id = start_frame_id

  if ( whole_line_num_1 > pre_base_block+blocks_num+end_base_block+pre_base+1 ):
    each = second_frame_id-start_frame_id
    time_step = (second_time-start_time)/each
    frames_num_2 = (end_frame_id-start_frame_id)/each+1
  else:
    frames_num_2 = 1
    time_step = 0.0
    each = 1

  if (frames_num_1 != frames_num_2):
    traj_tools.delete_duplicate(file_name, file_type)

  whole_line_num = len(open(file_name).readlines())
  frames_num = int((whole_line_num-pre_base)/(pre_base_block+blocks_num+end_base_block))

  #For groups, we will consider the connectivity.
  if return_group:
    if ( file_type == 'coord_xyz' or file_type == 'vel' or file_type == 'frc' ):

      element = []
      for i in range(blocks_num):
        line_i = linecache.getline(file_name, i+pre_base+pre_base_block+1)
        line_i_split = data_op.split_str(line_i, ' ')
        element.append(line_i_split[0])

      group_atom_1_id = []
      group_atoms_mass = []

      for i in range(len(group)):
        atom_id_i = atom_id[i]
        group_i_atom_1_id = []
        group_i_atoms_mass = []

        for j in range(len(group[i])):
          group_i_atoms_mass.append(atom.get_atom_mass(group[i][j])[1])
        group_atoms_mass.append(group_i_atoms_mass)

        for j in range(blocks_num-len(group[i])+1):
          if ( all( j+k+1 in atom_id_i and element[j+k] == group[i][k] for k in range(len(group[i]))) ):
            group_i_atom_1_id.append(j+1)
        group_atom_1_id.append(group_i_atom_1_id)

  linecache.clearcache()

  if ( file_type == 'coord_xyz' or file_type == 'coord_pdb' or file_type == 'vel' or file_type == 'frc' ):
    if return_group:
      return blocks_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, \
             end_frame_id, time_step, group_atom_1_id, group_atoms_mass
    else:
      return blocks_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step
  if ( file_type == 'ener' or file_type == 'mix_ener' ):
    return blocks_num, pre_base_block, end_base_block, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step
  if ( file_type == 'lagrange' ):
    return blocks_num, pre_base_block, end_base_block, pre_base, frames_num

if __name__ == '__main__':
  from CP2K_kit.tools import traj_info

  file_name = '/home/lujunbo/WORK/test/TEST_CP2K_KIT/TEST_ANALYZE/center/WATER_64H2O-pos-1.xyz'
  atom_id = [list(range(1,193))]
  groups = [['O','H','H']]
  blocks_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, \
  time_step, group_atom_1_id, group_atoms_mass = \
  traj_info.get_traj_info(file_name, 'coord_xyz', groups, atom_id, True)

  print (group_atom_1_id, group_atoms_mass)

