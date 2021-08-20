#! /usr/env/bin python

import os
import csv
import linecache
import subprocess
from collections import OrderedDict
from CP2K_kit.tools import list_dic_op
from CP2K_kit.tools import atom
from CP2K_kit.tools import call
from CP2K_kit.tools import get_cell
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import read_lmp
from CP2K_kit.tools import log_info
from CP2K_kit.deepff import lammps_run

hartree_to_ev = 27.2114
ang_to_bohr = 1.8897259886

def write_file(energy_cp2k, energy_lmp, frc_cp2k, frc_lmp, frc_x_cp2k, frc_x_lmp, \
               frc_y_cp2k, frc_y_lmp, frc_z_cp2k, frc_z_lmp, work_dir):

  '''
  write_file : write the energy and force data in a file

  Args :
    energy_cp2k : 1-d float list
      energy_cp2k contains cp2k energy values along with trajectory.
    energy_lmp : 1-d float list
      energy_lmp contains lammps energy values along with trajectory.
    frc_cp2k : 2-d float list, dim = Num of frames * (3*(Num of atoms))
      frc_cp2k contains cp2k force values along with trajectory.
    frc_lmp : 2-d float list, dim = Num of frames * (3*(Num of atoms))
      frc_lmp contains lammps force values along with trajectory.
    frc_x_cp2k : 2-d float list, dim = Num of frames * Num of atoms
      frc_x_cp2k contains cp2k force values of x part along with trajectory.
    frc_x_lmp : 2-d float list, dim = Num of frames * Num of atoms
      frc_x_lmp contains lammps force values of x part along with trajectory.
    frc_y_cp2k : 2-d float list, dim = Num of frames * Num of atoms
      frc_y_cp2k contains cp2k force values of y part along with trajectory.
    frc_y_lmp : 2-d float list, dim = Num of frames * Num of atoms
      frc_y_lmp contains lammps force values of y part along with trajectory.
    frc_z_cp2k : 2-d float list, dim = Num of frames * Num of atoms
      frc_z_cp2k contains cp2k force values of z part along with trajectory.
    frc_z_lmp : 2-d float list, dim = Num of frames * Num of atoms
      frc_z_lmp contains lammps force values of z part along with trajectory.
    work_dir : string
      working_dir is the working directory.
  Returns:
    none
  '''

  energy_file_name = ''.join((work_dir, '/energy.csv'))
  frc_file_name = ''.join((work_dir, '/force.csv'))
  frc_x_file_name = ''.join((work_dir, '/force_x.csv'))
  frc_y_file_name = ''.join((work_dir, '/force_y.csv'))
  frc_z_file_name = ''.join((work_dir, '/force_z.csv'))

  frames_num = len(energy_cp2k)
  atoms_num = len(frc_x_cp2k[0])

  with open(energy_file_name, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cp2k_energy (eV/atom)', 'lmp_energy (eV/atom)'])
    for i in range(frames_num):
      writer.writerow([energy_cp2k[i], energy_lmp[i]])

  with open(frc_file_name, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cp2k_force (eV/Angstrom)', 'lmp_force (eV/Angstrom)'])
    for i in range(frames_num):
      for j in range(atoms_num*3):
        writer.writerow([frc_cp2k[i][j], frc_lmp[i][j]])

  with open(frc_x_file_name, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cp2k_force (eV/Angstrom)', 'lmp_force (eV/Angstrom)'])
    for i in range(frames_num):
      for j in range(atoms_num):
        writer.writerow([frc_x_cp2k[i][j], frc_x_lmp[i][j]])

  with open(frc_y_file_name, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cp2k_force (eV/Angstrom)', 'lmp_force (eV/Angstrom)'])
    for i in range(frames_num):
      for j in range(atoms_num):
        writer.writerow([frc_y_cp2k[i][j], frc_y_lmp[i][j]])

  with open(frc_z_file_name, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['cp2k_force (eV/Angstrom)', 'lmp_force (eV/Angstrom)'])
    for i in range(frames_num):
      for j in range(atoms_num):
        writer.writerow([frc_z_cp2k[i][j], frc_z_lmp[i][j]])

def supervised_test(cp2k_pos_file, cp2k_cell_file, cp2k_frc_file, dpff_file, atom_label, work_dir):

  '''
  supervise_test : Do test for supervised learning

  Args :
    cp2k_pos_file : string
      cp2k_pos_file is the cp2k position trajectory file.
    cp2k_cell_file : string
      cp2k_cell_file is the cp2k cell trajectory file.
    cp2k_frc_file : string
      cp2k_frc_file is the cp2k force trajectory file.
    dpff_file : string
      dpff_file is the deepmd force field file.
    atom_label : dictionary
      atom_label is the atom label.
    work_dir : string
      work_dir is the working directory
  Returns:
    none
  '''

  atoms_num, base, pre_base, frames_num, each, start, end, time_step = traj_info.get_traj_info(cp2k_pos_file)

  tri_cell_tot = []
  atom_type_tot = []
  x_tot = []
  y_tot = []
  z_tot = []

  for i in range(frames_num):
    #Dump box
    line_i = linecache.getline(cp2k_cell_file, i+1+1)
    line_i_split = list_dic_op.str_split(line_i, ' ')
    cell_vec = [float(x) for x in line_i_split[2:11]]
    tri_cell_a, tri_cell_b, tri_cell_c = \
    get_cell.get_triclinic_cell(cell_vec[0:3], cell_vec[3:6], cell_vec[6:10])
    tri_cell_tot.append([str(tri_cell_a[0]), str(tri_cell_b[1]), str(tri_cell_c[2]), \
                         str(tri_cell_b[0]), str(tri_cell_c[0]), str(tri_cell_c[1])])

    #Dump coord
    atom_type_i = []
    x_i = []
    y_i = []
    z_i = []
    for j in range(atoms_num):
      line_ij = linecache.getline(cp2k_pos_file, (atoms_num+base)*i+j+pre_base+base+1)
      line_ij_split = list_dic_op.str_split(line_ij, ' ')
      line_ij_split[len(line_ij_split)-1] = line_ij_split[len(line_ij_split)-1].strip('\n')
      atom_type_i.append(list_dic_op.get_dic_keys(atom_label, line_ij_split[0])[0])
      x_i.append(line_ij_split[1])
      y_i.append(line_ij_split[2])
      z_i.append(line_ij_split[3])
    atom_type_tot.append(atom_type_i)
    x_tot.append(x_i)
    y_tot.append(y_i)
    z_tot.append(z_i)

  cmd = "cp %s %s" %(dpff_file, work_dir)
  call.call_simple_shell(work_dir, cmd)

  dpff_file_split = list_dic_op.str_split(dpff_file, '/')
  dpff_file_name = dpff_file_split[len(dpff_file_split)-1]

  for i in range(frames_num):
    frame_dir = ''.join((work_dir, '/frame_', str(i)))
    cmd = "mkdir %s" %(''.join(('frame_', str(i))))
    call.call_simple_shell(work_dir, cmd)
    lammps_run.gen_data_file(tri_cell_tot[i], atom_type_tot[i], x_tot[i], y_tot[i], z_tot[i], frame_dir, 'data.lmp')
    lmp_in_file_name = ''.join((frame_dir, '/in.lammps'))
    lmp_in_file = open(lmp_in_file_name, 'w')

    lmp_in_file.write('units           metal\n')
    lmp_in_file.write('boundary        p p p\n')
    lmp_in_file.write('atom_style      atomic\n')
    lmp_in_file.write('neighbor        1.0 bin\n')
    lmp_in_file.write('read_data       ./data.lmp\n')
    for key in atom_label.keys():
      mass = atom.get_atom_mass(atom_label[key])[1]
      lmp_in_file.write('mass            %d %f\n' %(key, mass))

    lmp_in_file.write('pair_style      deepmd %s\n' %(''.join((work_dir, '/', dpff_file_name))))
    lmp_in_file.write('pair_coeff\n')

    lmp_in_file.write('thermo_style    custom step temp pe ke etotal\n')
    lmp_in_file.write('thermo          1\n')
    lmp_in_file.write('dump            1 all custom 1 atom.dump id type x y z fx fy fz\n')
    lmp_in_file.write('run 0\n')

    lmp_in_file.close()

    cmd = 'lmp < ./in.lammps > out'
    subprocess.run(cmd, shell=True, cwd=frame_dir)

  energy_cp2k = []
  energy_lmp = []
  frc_cp2k = []
  frc_lmp = []
  frc_x_cp2k = []
  frc_x_lmp = []
  frc_y_cp2k = []
  frc_y_lmp = []
  frc_z_cp2k = []
  frc_z_lmp = []

  for i in range(frames_num):
    frame_dir = ''.join((work_dir, '/frame_', str(i)))
    lmp_log_file = ''.join((frame_dir, '/log.lammps'))
    lmp_traj_file = ''.join((frame_dir, '/atom.dump'))
    cmd_a = "grep -n %s %s" % ("'Step'", 'log.lammps')
    a = call.call_returns_shell(frame_dir, cmd_a)
    a_int = int(a[0].split(':')[0])
    line_log_i = linecache.getline(lmp_log_file, a_int+1)
    line_log_i_split = list_dic_op.str_split(line_log_i, ' ')
    energy_lmp.append(float(line_log_i_split[2])/atoms_num)
    atom_id_i = []
    frc_lmp_i = []
    frc_x_lmp_i = []
    frc_y_lmp_i = []
    frc_z_lmp_i = []
    for j in range(atoms_num):
      line_traj_ij = linecache.getline(lmp_traj_file, j+1+9)
      line_traj_ij_split = list_dic_op.str_split(line_traj_ij, ' ')
      atom_id_i.append(int(line_traj_ij_split[0]))
      frc_x_lmp_i.append(float(line_traj_ij_split[5]))
      frc_y_lmp_i.append(float(line_traj_ij_split[6]))
      frc_z_lmp_i.append(float(line_traj_ij_split[7].strip('\n')))
      frc_lmp_i.append([float(line_traj_ij_split[5]), float(line_traj_ij_split[6]), float(line_traj_ij_split[7].strip('\n'))])

    atom_id_i_asc, asc_index = list_dic_op.list_order(atom_id_i, 'ascend', True)
    frc_x_lmp_i_asc = list_dic_op.order_list(frc_x_lmp_i, asc_index)
    frc_y_lmp_i_asc = list_dic_op.order_list(frc_y_lmp_i, asc_index)
    frc_z_lmp_i_asc = list_dic_op.order_list(frc_z_lmp_i, asc_index)
    frc_lmp_i_asc = list_dic_op.order_list(frc_lmp_i, asc_index)

    frc_lmp.append(list_dic_op.list_reshape(frc_lmp_i_asc))
    frc_x_lmp.append(frc_x_lmp_i_asc)
    frc_y_lmp.append(frc_y_lmp_i_asc)
    frc_z_lmp.append(frc_z_lmp_i_asc)

    line_cp2k = linecache.getline(cp2k_frc_file, i*(atoms_num+base)+2)
    line_cp2k_split = list_dic_op.str_split(line_cp2k, ' ')
    energy_cp2k.append(float(line_cp2k_split[len(line_cp2k_split)-1].strip('\n'))*hartree_to_ev/atoms_num)

    frc_cp2k_i = []
    frc_x_cp2k_i = []
    frc_y_cp2k_i = []
    frc_z_cp2k_i = []
    for j in range(atoms_num):
      line_cp2k_ij = linecache.getline(cp2k_frc_file, i*(atoms_num+base)+j+1+pre_base+base)
      line_cp2k_ij_split = list_dic_op.str_split(line_cp2k_ij, ' ')
      fx = float(line_cp2k_ij_split[1])*hartree_to_ev*ang_to_bohr
      fy = float(line_cp2k_ij_split[2])*hartree_to_ev*ang_to_bohr
      fz = float(line_cp2k_ij_split[3].strip('\n'))*hartree_to_ev*ang_to_bohr
      frc_x_cp2k_i.append(fx)
      frc_y_cp2k_i.append(fy)
      frc_z_cp2k_i.append(fz)
      frc_cp2k_i.append(fx)
      frc_cp2k_i.append(fy)
      frc_cp2k_i.append(fz)

    frc_cp2k.append(frc_cp2k_i)
    frc_x_cp2k.append(frc_x_cp2k_i)
    frc_y_cp2k.append(frc_y_cp2k_i)
    frc_z_cp2k.append(frc_z_cp2k_i)

  write_file(energy_cp2k, energy_lmp, frc_cp2k, frc_lmp, frc_x_cp2k, frc_x_lmp, \
             frc_y_cp2k, frc_y_lmp, frc_z_cp2k, frc_z_lmp, work_dir)

def active_learning_test(lmp_traj_file, lmp_log_file, cp2k_inp_file, atom_label, work_dir):

  '''
  active_learning_test : Do test for active learning

  Args :
    lmp_traj_file : string
      lmp_traj_file is the lammps trajectory file.
    lmp_log_file : string
      lmp_log_file is the lammps output file.
    cp2k_inp_file : string
      cp2k_inp_file is the cp2k input file.
    atom_label : dictionary
      atom_label is the atom label.
    work_dir : string
      work_dir is the workding directory.
  Returns :
    none
  '''

  energy_cp2k = []
  frc_cp2k = []
  frc_x_cp2k = []
  frc_y_cp2k = []
  frc_z_cp2k = []
  energy_lmp = []
  frc_lmp = []
  frc_x_lmp = []
  frc_y_lmp = []
  frc_z_lmp = []

  atoms_num, frames_num, start_id, end_id, each = read_lmp.lmp_traj_info(lmp_traj_file, lmp_log_file)

  atoms, energy, coord, vel, frc, cell = \
  read_lmp.read_lmp_log_traj(lmp_traj_file, lmp_log_file, atom_label, [], True, True, False, True, True)

  cp2k_inp_split = list_dic_op.str_split(os.path.abspath(cp2k_inp_file), '/')
  cp2k_inp_file = cp2k_inp_split[len(cp2k_inp_split)-1]

  for i in range(frames_num):
    frame_dir = ''.join((work_dir, '/frame_', str(i)))
    cmd = "mkdir %s" %(''.join(('frame_', str(i))))
    call.call_simple_shell(work_dir, cmd)

    cmd = "cp %s %s" %(cp2k_inp_file, ''.join(('frame_', str(i))))
    call.call_simple_shell(work_dir, cmd)

    cmd = "cp %s %s" %(cp2k_inp_file ,frame_dir)
    call.call_simple_shell(work_dir, cmd)

    cell_file_name = ''.join((frame_dir, '/box.inc'))
    cell_file = open(cell_file_name, 'w')
    cell_file.write('%s    %f    %f    %f\n' %('A', cell[i][0][0], cell[i][0][1], cell[i][0][2]))
    cell_file.write('%s    %f    %f    %f\n' %('B', cell[i][1][0], cell[i][1][1], cell[i][1][2]))
    cell_file.write('%s    %f    %f    %f\n' %('C', cell[i][2][0], cell[i][2][1], cell[i][2][2]))
    cell_file.close()
    coord_file_name = ''.join((frame_dir, '/coord.inc'))
    coord_file = open(coord_file_name, 'w')
    for j in range(atoms_num):
      coord_file.write('%s    %f    %f    %f\n' %(atoms[i][j], coord[i][j][0], coord[i][j][1], coord[i][j][2]))
    coord_file.close()

    cmd = "cp2k.popt %s 1> cp2k.out 2> cp2k.err" %(cp2k_inp_file)
    subprocess.run(cmd, shell=True, cwd=frame_dir)

  for i in range(frames_num):
    frame_dir = ''.join((work_dir, '/frame_', str(i)))
    frc_cp2k_i = []
    frc_x_cp2k_i = []
    frc_y_cp2k_i = []
    frc_z_cp2k_i = []
    frc_lmp_i = []
    frc_x_lmp_i = []
    frc_y_lmp_i = []
    frc_z_lmp_i = []

    energy_lmp.append(energy[i]/atoms_num)
    cmd = "grep %s %s" % ("'ENERGY| Total FORCE_EVAL'", 'cp2k.out')
    energy_parse = call.call_returns_shell(frame_dir, cmd)
    energy_cp2k_i = float(energy_parse[0].split(':')[1].strip('\n'))
    energy_cp2k_i = energy_cp2k_i*hartree_to_ev/atoms_num
    energy_cp2k.append(energy_cp2k_i)

    frc_file_i = ''.join((frame_dir, '/cp2k-1_0.xyz'))
    for j in range(atoms_num):
      line_ij = linecache.getline(frc_file_i, j+4+1)
      line_ij_split = list_dic_op.str_split(line_ij, ' ')
      fx = float(line_ij_split[3])*hartree_to_ev*ang_to_bohr
      fy = float(line_ij_split[4])*hartree_to_ev*ang_to_bohr
      fz = float(line_ij_split[5].strip('\n'))*hartree_to_ev*ang_to_bohr
      frc_cp2k_i.append(fx)
      frc_cp2k_i.append(fy)
      frc_cp2k_i.append(fz)
      frc_x_cp2k_i.append(fx)
      frc_y_cp2k_i.append(fy)
      frc_z_cp2k_i.append(fz)

      frc_lmp_i.append(frc[i][j][0])
      frc_lmp_i.append(frc[i][j][1])
      frc_lmp_i.append(frc[i][j][2])
      frc_x_lmp_i.append(frc[i][j][0])
      frc_y_lmp_i.append(frc[i][j][1])
      frc_z_lmp_i.append(frc[i][j][2])

    frc_cp2k.append(frc_cp2k_i)
    frc_x_cp2k.append(frc_x_cp2k_i)
    frc_y_cp2k.append(frc_y_cp2k_i)
    frc_z_cp2k.append(frc_z_cp2k_i)

    frc_lmp.append(frc_lmp_i)
    frc_x_lmp.append(frc_x_lmp_i)
    frc_y_lmp.append(frc_y_lmp_i)
    frc_z_lmp.append(frc_z_lmp_i)

  write_file(energy_cp2k, energy_lmp, frc_cp2k, frc_lmp, frc_x_cp2k, frc_x_lmp, \
             frc_y_cp2k, frc_y_lmp, frc_z_cp2k, frc_z_lmp, work_dir)

def dp_test_run(dp_test_param, work_dir):

  '''
  dp_test_run : kernel function to run dp_test

  Args :
    dp_test_param : dictionary
      dp_test_param contains keywords used in deep potential test.
    work_dir : string
      work_dir is the working directory.
  Returns :
    none
  '''

  if ( 'learn_type' in dp_test_param.keys() ):
    learn_type = dp_test_param['learn_type']
  else:
    log_info.log_error('No learning type found, please set analyze/dp_test/learn_type')
    exit()

  if ( learn_type == 'supervised' ):
    if ( 'cp2k_frc_file' in dp_test_param.keys() ):
      cp2k_frc_file = dp_test_param['cp2k_frc_file']
    else:
      log_info.log_error('No cp2k force trajectory found, please set analyze/dp_test/cp2k_frc_file')
      exit()
    if ( 'cp2k_pos_file' in dp_test_param.keys() ):
      cp2k_pos_file = dp_test_param['cp2k_pos_file']
    else:
      log_info.log_error('No cp2k position trajectory found, please set analyze/dp_test/cp2k_pos_file')
      exit()
    if ( 'cp2k_cell_file' in dp_test_param.keys() ):
      cp2k_cell_file = dp_test_param['cp2k_cell_file']
    else:
      log_info.log_error('No cp2k cell trajectory found, please set analyze/dp_test/cp2k_cell_file')
      exit()
    if ( 'dpff_file' in dp_test_param.keys() ):
      dpff_file = dp_test_param['dpff_file']
    else:
      log_info.log_error('No deepmd-kit force field file found, please set analyze/dp_test/dpff_file')
      exit()
    if ( 'atom_label' in dp_test_param.keys() ):
      atom_label = dp_test_param['atom_label']
      atom_label_dic = OrderedDict()
      for i in range (len(atom_label)):
        lable_split = list_dic_op.str_split(atom_label[i], ':')
        atom_label_dic[int(lable_split[0])] = lable_split[1]
    else:
      log_info.log_error('No atom lable found, please set analyze/dp_test/atom_lable')
      exit()

    supervised_test(cp2k_pos_file, cp2k_cell_file, cp2k_frc_file, dpff_file, atom_label_dic, work_dir)

  elif ( learn_type == 'active_learning' ):
    if ( 'lmp_traj_file' in dp_test_param.keys() ):
      lmp_traj_file = dp_test_param['lmp_traj_file']
    else:
      log_info.log_error('No lammps trajectory found, please set analyze/dp_test/lmp_traj_file')
      exit()
    if ( 'lmp_log_file' in dp_test_param.keys() ):
      lmp_log_file = dp_test_param['lmp_log_file']
    else:
      log_info.log_error('No lammps output file found, please set analyze/dp_test/lmp_log_file')
      exit()
    if ( 'cp2k_inp_file' in dp_test_param.keys() ):
      cp2k_inp_file = dp_test_param['cp2k_inp_file']
    else:
      log_info.log_error('No cp2k input file found, please set analyze/dp_test/cp2k_inp_file')
      exit()
    if ( 'atom_label' in dp_test_param.keys() ):
      atom_label = dp_test_param['atom_label']
      atom_label_dic = OrderedDict()
      for i in range (len(atom_label)):
        lable_split = list_dic_op.str_split(atom_label[i], ':')
        atom_label_dic[int(lable_split[0])] = lable_split[1]
    else:
      log_info.log_error('No atom label found, please set analyze/dp_test/atom_lable')
      exit()
    active_learning_test(lmp_traj_file, lmp_log_file, cp2k_inp_file, atom_label_dic, work_dir)

  else:
    log_info.log_error('Could not recognize learn_type, please check')
    exit()
