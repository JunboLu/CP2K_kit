#! /usr/env/bin python

import os
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import get_cell
from CP2K_kit.tools import data_op

def lmp_traj_info(lmp_traj_file, lmp_log_file):

  '''
  lmp_traj_info: Dump information from lammps trajectory file and output file.

  Args:
    lmp_traj_file: string
      lmp_traj_file is the lammps trajectory file.
    lmp_log_file: string
      lmp_log_file is the lammps output file.
  Returns:
    atoms_num: int
      atoms_num is the number of atoms for one frame in lammps trajectory.
    frames_num: int
      frames_num is the number of frames in lammps trajectory.
    start_id: int
      start_id is the starting index.
    end_id: int
      end_id is the endding index.
    each: int
      each is the increment.
  '''

  line = linecache.getline(lmp_traj_file, 4)
  atoms_num = int(line.strip('\n'))
  traj_frames_num = int(len(open(lmp_traj_file).readlines())/(atoms_num+9))

  cmd_a = "grep -n %s %s" % ("'Step'", os.path.abspath(lmp_log_file))
  a = call.call_returns_shell(os.getcwd(), cmd_a)
  a_int = int(a[0].split(':')[0])

  line = linecache.getline(lmp_log_file, a_int)
  line_split = data_op.split_str(line, ' ', '\n')
  log_id_num = len(line_split)
  step_id = line_split.index('Step')

  line = linecache.getline(lmp_log_file, a_int+1)
  line_split = data_op.split_str(line, ' ')
  start_id = int(line_split[step_id])

  line = linecache.getline(lmp_log_file, a_int+2)
  line_split = data_op.split_str(line, ' ', '\n')
  if ( len(line_split) == log_id_num ):
    each = int(line_split[step_id])-start_id
  else:
    each = 0

  line_num = a_int+traj_frames_num
  while True:
    line = linecache.getline(lmp_log_file, line_num)
    line_split = data_op.split_str(line, ' ', '\n')
    if ( len(line_split) != log_id_num ):
      line_num = line_num-1
    else:
      break

  line = linecache.getline(lmp_log_file, line_num)
  line_split = data_op.split_str(line, ' ', '\n')
  end_id = int(line_split[step_id])
  if ( each != 0 ):
    frames_num = int((end_id-start_id)/each)+1
  else:
    frames_num = 1

  linecache.clearcache()

  return atoms_num, frames_num, start_id, end_id, each

def read_lmp_log_traj(lmp_traj_file, lmp_log_file, atom_label={}, frames=[], ene_return=False, \
                      coord_return=False, vel_return=False, frc_return=False, cell_return=False):

  '''
  read_lmp_log_traj: read lammps trajectory and output files.

  Args:
    lmp_traj_file: string
      lmp_traj_file is the lammps trajectory file.
    lmp_log_file: string
      lmp_log_file is the lammps output file.
    atom_label: dictionary
      atom_label is the atomic label.
    frames: 1-d int list
      frames is the choosed frame_id
    ene_return: bool
      ene is whether we need to return energy.
    coord_return: bool
      ene is whether we need to return coordinates.
    vel_return: bool
      vel is whether we need to return velocities.
    frc_return: bool
      ene is whether we need to return forces.
    cell_return: bool
      ene is whether we need to return cells.
  Returns:
    atoms: 2-d string list, dim = (num of frames)*(num of atoms)
      atoms is the atom name along with the trajectory.
    energy: 2-d float list, dim = (num of frames)*(num of atoms)
      energy is the energy along with the trajectory.
    coord: 3-d float list, dim = (num of frames)*(num of atoms)*3
      coord is the coordinates along with the trajectory.
    vel: 3-d float list, dim = (num of frames)*(num of atoms)*3
      vel is the velocity along with the trajectory.
    frc: 3-d float list, dim = (num of frames)*(num of atoms)*3
      frc is the force along with the trajectory.
    cell: 3-d float list, dim = (num of frames)*3*3
      cell is the cell vector along with the trajectory.
  '''

  line = linecache.getline(lmp_traj_file, 9)
  line_split = data_op.split_str(line, ' ', '\n')
  line_split[len(line_split)-1] = line_split[len(line_split)-1]

  if ( 'id' in line_split ):
    atom_id_id = line_split.index('id')-2
  else:
    print ('Could not find id in dump item, please check!')
    exit()
  if ( 'type' in line_split ):
    atom_type_id = line_split.index('type')-2
  else:
    print ('Could not find atom type id in dump item, please check!')
    exit()
  if coord_return:
    if ( 'x' in line_split and 'y' in line_split and 'z' in line_split ):
      x_id = line_split.index('x')-2
      y_id = line_split.index('y')-2
      z_id = line_split.index('z')-2
    else:
      print ('You want the coordinates, but could not find xyz id in dump item, please check!')
      exit()
  if vel_return:
    if ( 'vx' in line_split and 'vy' in line_split and 'vz' in line_split):
      vx_id = line_split.index('vx')-2
      vy_id = line_split.index('vy')-2
      vz_id = line_split.index('vz')-2
    else:
      print ('You want the velocity, but could not find vxvyvz id in dump item, please check!')
      exit()
  if frc_return:
    if ( 'fx' in line_split and 'fy' in line_split and 'fz' in line_split):
      fx_id = line_split.index('fx')-2
      fy_id = line_split.index('fy')-2
      fz_id = line_split.index('fz')-2
    else:
      print ('You want the force, but could not find fxfyfz id in dump item, please check!')
      exit()

  atoms = []
  energy = []
  coord = []
  vel = []
  frc = []
  cell = []

  cmd_a = "grep -n %s %s" % ("'Step'", lmp_log_file)
  a = call.call_returns_shell(os.getcwd(), cmd_a)
  a_int = int(a[0].split(':')[0])
  line = linecache.getline(lmp_log_file, a_int)
  line_split = data_op.split_str(line, ' ', '\n')
  line_split[len(line_split)-1] = line_split[len(line_split)-1]

  if ene_return:
    if ( 'PotEng' in line_split ):
      ene_id = line_split.index('PotEng')
    else:
      print ('You want the potential energy, but could not find PotEng id in thermo item, please check!')
      exit()

  atoms_num, frames_num, start_id, end_id, each = lmp_traj_info(lmp_traj_file, lmp_log_file)

  if ( frames == [] ):
    for i in range(frames_num):
      frames.append(start_id+each*i)

  for i in frames:
    if ene_return:
      line_log_i = linecache.getline(lmp_log_file, int((i-start_id)/each)+a_int+1)
      line_log_i_split = data_op.split_str(line_log_i, ' ')
      energy.append(float(line_log_i_split[ene_id]))

    if cell_return:
      line_1 = linecache.getline(lmp_traj_file, (atoms_num+9)*int((i-start_id)/each)+6)
      line_1_split = data_op.split_str(line_1, ' ', '\n')
      line_2 = linecache.getline(lmp_traj_file, (atoms_num+9)*int((i-start_id)/each)+7)
      line_2_split = data_op.split_str(line_2, ' ', '\n')
      line_3 = linecache.getline(lmp_traj_file, (atoms_num+9)*int((i-start_id)/each)+8)
      line_3_split = data_op.split_str(line_3, ' ', '\n')
      Lx = float(line_1_split[1]) - float(line_1_split[0])
      Ly = float(line_2_split[1]) - float(line_2_split[0])
      Lz = float(line_3_split[1]) - float(line_3_split[0])
      xy = float(line_1_split[2])
      xz = float(line_2_split[2])
      yz = float(line_3_split[2])
      a_vec, b_vec, c_vec = get_cell.get_triclinic_cell_six([Lx, Ly, Lz, xy, xz, yz])
      cell.append([a_vec, b_vec, c_vec])

    atom_type_i = []
    atom_id_i = []
    frc_i = []
    coord_i = []
    vel_i = []

    for j in range(atoms_num):
      line_ij = linecache.getline(lmp_traj_file, (atoms_num+9)*int((i-start_id)/each)+j+1+9)
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      line_ij_split[len(line_ij_split)-1] = line_ij_split[len(line_ij_split)-1]
      atom_type = int(line_ij_split[atom_type_id])
      atom_id = int(line_ij_split[atom_id_id])
      atom_id_i.append(atom_id)
      if ( atom_label != {} ):
        atom_type_i.append(atom_label[atom_type])
      if coord_return:
        x = float(line_ij_split[x_id])
        y = float(line_ij_split[y_id])
        z = float(line_ij_split[z_id])
        coord_i.append([x,y,z])
      if vel_return:
        vx = float(line_ij_split[vx_id])
        vy = float(line_ij_split[vy_id])
        vz = float(line_ij_split[vz_id])
        vel_i.append([vx,vy,vz])
      if frc_return:
        fx = float(line_ij_split[fx_id])
        fy = float(line_ij_split[fy_id])
        fz = float(line_ij_split[fz_id])
        frc_i.append([fx,fy,fz])

    atom_id_i_asc, asc_index = data_op.get_list_order(atom_id_i, 'ascend', True)
    if ( atom_label != {} ):
      atom_type_i_asc = data_op.reorder_list(atom_type_i, asc_index)
      atoms.append(atom_type_i_asc)
    if coord_return:
      coord_i_asc = data_op.reorder_list(coord_i, asc_index)
      coord.append(coord_i_asc)
    if vel_return:
      vel_i_asc = data_op.reorder_list(vel_i, asc_index)
      vel.append(vel_i_asc)
    if frc_return:
      frc_i_asc = data_op.reorder_list(frc_i, asc_index)
      frc.append(frc_i_asc)

  linecache.clearcache()

  return atoms, energy, coord, vel, frc, cell
