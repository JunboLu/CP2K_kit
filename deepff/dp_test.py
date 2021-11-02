#! /usr/env/bin python

import os
import linecache
import subprocess
from collections import OrderedDict
from CP2K_kit.tools import data_op
from CP2K_kit.tools import atom
from CP2K_kit.tools import call
from CP2K_kit.tools import get_cell
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import read_lmp
from CP2K_kit.tools import log_info
from CP2K_kit.deepff import lammps_run
from CP2K_kit.deepff import model_devi
from CP2K_kit.deepff import check_deepff
from CP2K_kit.deepff import write_data
from CP2K_kit.deepff import process

hartree_to_ev = 2.72113838565563E+01
ang_to_bohr = 1.0/5.29177208590000E-01

def supervised_test(cp2k_pos_file, cp2k_cell_file, cp2k_frc_file, lmp_exe, lmp_mpi_num, lmp_omp_num, dpff_file, atom_label, work_dir):

  '''
  supervise_test: Do test for supervised learning

  Args :
    cp2k_pos_file: string
      cp2k_pos_file is the cp2k position trajectory file.
    cp2k_cell_file: string
      cp2k_cell_file is the cp2k cell trajectory file.
    cp2k_frc_file: string
      cp2k_frc_file is the cp2k force trajectory file.
    lmp_exe: string
      lmp_exe is the lammps executable file.
    lmp_mpi_num: int
      lmp_mpi_num is the lammps mpi number.
    lmp_omp_num: int
      lmp_omp_num is the lammps openmp number.
    dpff_file: string
      dpff_file is the deepmd force field file.
    atom_label: dictionary
      atom_label is the atom label.
    work_dir: string
      work_dir is the working directory
  Returns:
    none
  '''

  atoms_num, pre_base_block, end_base_block, pre_base, frames_num, each, start, end, time_step = \
  traj_info.get_traj_info(cp2k_pos_file, 'coord_xyz')

  tri_cell_tot = []
  atom_type_tot = []
  x_tot = []
  y_tot = []
  z_tot = []

  for i in range(frames_num):
    #Dump box
    line_i = linecache.getline(cp2k_cell_file, i+1+1)
    line_i_split = data_op.split_str(line_i, ' ')
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
      line_ij = linecache.getline(cp2k_pos_file, (pre_base_block+atoms_num+end_base_block)*i+j+pre_base+pre_base_block+1)
      line_ij_split = data_op.split_str(line_ij, ' ', '\n')
      atom_type_i.append(data_op.get_dic_keys(atom_label, line_ij_split[0])[0])
      x_i.append(line_ij_split[1])
      y_i.append(line_ij_split[2])
      z_i.append(line_ij_split[3])
    atom_type_tot.append(atom_type_i)
    x_tot.append(x_i)
    y_tot.append(y_i)
    z_tot.append(z_i)

  linecache.clearcache()

  cmd = "cp %s %s" %(dpff_file, work_dir)
  call.call_simple_shell(work_dir, cmd)

  dpff_file_split = data_op.split_str(dpff_file, '/')
  dpff_file_name = dpff_file_split[len(dpff_file_split)-1]

  print ('Run lammps jobs for %s system' %(frames_num), flush=True)

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

    cmd = 'export OMP_NUM_THREADS=%d && mpirun -np %d %s < ./in.lammps > out' %(lmp_omp_num, lmp_mpi_num, lmp_exe)
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
    line_log_i_split = data_op.split_str(line_log_i, ' ')
    energy_lmp.append(float(line_log_i_split[2])/atoms_num)
    atom_id_i = []
    frc_lmp_i = []
    frc_x_lmp_i = []
    frc_y_lmp_i = []
    frc_z_lmp_i = []
    for j in range(atoms_num):
      line_traj_ij = linecache.getline(lmp_traj_file, j+1+9)
      line_traj_ij_split = data_op.split_str(line_traj_ij, ' ', '\n')
      atom_id_i.append(int(line_traj_ij_split[0]))
      frc_x_lmp_i.append(float(line_traj_ij_split[5]))
      frc_y_lmp_i.append(float(line_traj_ij_split[6]))
      frc_z_lmp_i.append(float(line_traj_ij_split[7]))
      frc_lmp_i.append([float(line_traj_ij_split[5]), float(line_traj_ij_split[6]), float(line_traj_ij_split[7])])

    atom_id_i_asc, asc_index = data_op.get_list_order(atom_id_i, 'ascend', True)
    frc_x_lmp_i_asc = data_op.order_list(frc_x_lmp_i, asc_index)
    frc_y_lmp_i_asc = data_op.order_list(frc_y_lmp_i, asc_index)
    frc_z_lmp_i_asc = data_op.order_list(frc_z_lmp_i, asc_index)
    frc_lmp_i_asc = data_op.order_list(frc_lmp_i, asc_index)

    frc_lmp.append(data_op.list_reshape(frc_lmp_i_asc))
    frc_x_lmp.append(frc_x_lmp_i_asc)
    frc_y_lmp.append(frc_y_lmp_i_asc)
    frc_z_lmp.append(frc_z_lmp_i_asc)

    line_cp2k = linecache.getline(cp2k_frc_file, i*(pre_base_block+atoms_num+end_base_block)+2+pre_base)
    line_cp2k_split = data_op.split_str(line_cp2k, ' ', '\n')
    energy_cp2k.append(float(line_cp2k_split[len(line_cp2k_split)-1])*hartree_to_ev/atoms_num)

    frc_cp2k_i = []
    frc_x_cp2k_i = []
    frc_y_cp2k_i = []
    frc_z_cp2k_i = []
    for j in range(atoms_num):
      line_cp2k_ij = linecache.getline(cp2k_frc_file, i*(pre_base_block+atoms_num+end_base_block)+j+1+pre_base+pre_base_block)
      line_cp2k_ij_split = data_op.split_str(line_cp2k_ij, ' ', '\n')
      fx = float(line_cp2k_ij_split[1])*hartree_to_ev*ang_to_bohr
      fy = float(line_cp2k_ij_split[2])*hartree_to_ev*ang_to_bohr
      fz = float(line_cp2k_ij_split[3])*hartree_to_ev*ang_to_bohr
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

  linecache.clearcache()

  write_data.write_file(energy_cp2k, energy_lmp, frc_cp2k, frc_lmp, frc_x_cp2k, frc_x_lmp, \
                        frc_y_cp2k, frc_y_lmp, frc_z_cp2k, frc_z_lmp, work_dir)

def active_learning_test(work_dir, iter_id, atoms_type_multi_sys, use_mtd_tot, force_conv, energy_conv):

  '''
  active_learning_test: Do test for active learning

  Args:
    work_dir: string
      work_dir is workding directory.
    iter_id: int
      iter_id is current iteration number.
    atom_type_dic_tot: dictionary
      Example: {0:{'O':1,'H':2}, 1:{'O':1,'H':2}}, The keys stands for system.
    use_mtd_tot: 1-d bool list
      use_mtd_tot is whether using metadynamics for different systems.
    force_conv: float
      force_conv is the maximum force convergence.
    energy_conv: float
      energy_conv is the maximum energy convergence.
  Returns:
    choosed_index:
  '''

  struct_index = OrderedDict()
  success_frames = []
  tot_frames = []

  lmp_dir = ''.join((work_dir, '/iter_', str(iter_id), '/02.lammps_calc'))
  sys_num = len(atoms_type_multi_sys)
  for i in range(sys_num):
    struct_index_i = OrderedDict()
    success_frames_i = []
    tot_frames_i = []
    lmp_sys_dir = ''.join((lmp_dir, '/sys_', str(i)))
    task_num = process.get_task_num(lmp_sys_dir)
    for j in range(task_num):
      lmp_sys_task_dir = ''.join((lmp_sys_dir, '/task_', str(j)))
      lmp_log_file = ''.join((lmp_sys_task_dir, '/lammps.out'))
      lmp_traj_file = ''.join((lmp_sys_task_dir, '/atom.dump'))
      atoms_num, frames_num, start_id, end_id, each = read_lmp.lmp_traj_info(lmp_traj_file, lmp_log_file)
      atoms, energy_lmp, coord_lmp, vel_lmp, frc_lmp, cell_lmp = \
      read_lmp.read_lmp_log_traj(lmp_traj_file, lmp_log_file, atoms_type_multi_sys[i], [], True, True, False, True, True)
      tot_frames_i.append(frames_num)
      if ( use_mtd_tot[i] == 'mtd' ):
        frc_lmp = []
        model_dir = ''.join((lmp_sys_task_dir, '/model_0'))
        for k in range(frames_num):
          model_traj_dir = ''.join((model_dir, '/traj_', str(k)))
          id_k = []
          frc_k = []
          dump_file = ''.join((model_traj_dir, '/atom.dump'))
          for l in range(atoms_num):
            line_kl = linecache.getline(dump_file, l+9+1)
            line_kl_split = data_op.split_str(line_kl, ' ', '\n')
            id_k.append(int(line_lm_split[0]))
            frc_k.append([float(line_kl_split[5]), float(line_kl_split[6]), float(line_kl_split[7])])
          linecache.clearcache()
          id_k_asc, asc_index = data_op.get_list_order(id_k, 'ascend', True)
          frc_k_asc = data_op.reorder_list(frc_k, asc_index)
          frc_lmp.append(frc_k_asc)

      energy_cp2k_final = []
      energy_lmp_final = []
      frc_cp2k_final = []
      frc_lmp_final = []
      index_final = []
      cp2k_sys_task_dir = ''.join((work_dir, '/iter_', str(iter_id), '/03.cp2k_calc/sys_', str(i), '/task_', str(j)))
      for k in range(frames_num):
        cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(k)))
        frc_cp2k_k = []
        frc_lmp_k = []
        cmd = "grep %s %s" % ("'ENERGY| Total FORCE_EVAL'", 'cp2k.out')
        energy_parse = call.call_returns_shell(cp2k_sys_task_traj_dir, cmd)
        frc_file_k = ''.join((cp2k_sys_task_traj_dir, '/cp2k-1_0.xyz'))
        if ( len(energy_parse) != 0 and os.path.exists(frc_file_k) ):
          index_final.append(k)
          energy_lmp_final.append(energy_lmp[k]/atoms_num)
          energy_cp2k_k = float(energy_parse[0].split(':')[1].strip('\n'))
          energy_cp2k_k = energy_cp2k_k*hartree_to_ev/atoms_num
          energy_cp2k_final.append(energy_cp2k_k)

          for l in range(atoms_num):
            line_kl = linecache.getline(frc_file_k, l+4+1)
            line_kl_split = data_op.split_str(line_kl, ' ', '\n')
            fx = float(line_kl_split[3])*hartree_to_ev*ang_to_bohr
            fy = float(line_kl_split[4])*hartree_to_ev*ang_to_bohr
            fz = float(line_kl_split[5])*hartree_to_ev*ang_to_bohr
            frc_cp2k_k.append([fx, fy, fz])
            frc_lmp_k.append([frc_lmp[k][l][0], frc_lmp[k][l][1], frc_lmp[k][l][0]])
          linecache.clearcache()

          frc_cp2k_final.append(frc_cp2k_k)
          frc_lmp_final.append(frc_lmp_k)

      write_data.write_file(energy_cp2k_final, energy_lmp_final, frc_cp2k_final, frc_lmp_final, cp2k_sys_task_dir)

      choosed_index = []
      success_frames_ij = 0
      for k in range(len(energy_cp2k_final)):
        #Be careful, we use index k for force_devi.
        force_devi = model_devi.calc_force_devi(frc_cp2k_final[k], frc_lmp_final[k])
        #Be careful, we use index_final[k] for atoms_type_dist.
        atoms_type_dist = model_devi.calc_dist(atoms, coord_lmp[index_final[k]], cell_lmp[index_final[k]][0], \
                                          cell_lmp[index_final[k]][1], cell_lmp[index_final[k]][2])
        dist = []
        atom_type_pair_tot = []
        for key in atoms_type_dist:
          dist.append(min(atoms_type_dist[key]))
          atom_type_pair_tot.append(key)
        min_dist = min(dist)
        min_dist_index = dist.index(min_dist)
        atom_type_pair = atom_type_pair_tot[min_dist_index]

        if ( max(force_devi) < force_conv and abs(energy_cp2k_final[i]-energy_lmp_final[i]) < energy_conv):
          success_frames_ij = success_frames_ij+1
        else:
          atom_cov_radii_plus = atom.get_atom_cov_radius(atom_type_pair[0]) + \
                                atom.get_atom_cov_radius(atom_type_pair[1])
          if ( min_dist > atom_cov_radii_plus*0.7 ):
            choosed_index.append(index_final[k])
      success_frames_i.append(success_frames_ij)
      struct_index_i[j] = choosed_index
    tot_frames.append(tot_frames_i)
    success_frames.append(success_frames_i)
    struct_index[i] = struct_index_i

  success_ratio_sys = []
  tot_frames_sys = []
  success_frames_sys = []
  for i in range(len(tot_frames)):
    success_ratio_sys.append(float(float(sum(success_frames[i]))/float(sum(tot_frames[i]))))
    tot_frames_sys.append(sum(tot_frames[i]))
    success_frames_sys.append(sum(success_frames[i]))

  success_ratio = float(float(sum(success_frames_sys))/float(sum(tot_frames_sys)))

  return struct_index, success_ratio_sys, success_ratio

def dp_test_run(dp_test_param, work_dir):

  '''
  dp_test_run: kernel function to run dp_test

  Args:
    dp_test_param: dictionary
      dp_test_param contains keywords used in deep potential test.
    work_dir: string
      work_dir is the working directory.
  Returns:
    none
  '''

  dp_test_param = check_deepff.check_dp_test(dp_test_param)

  learn_type = dp_test_param['learn_type']
  atom_label = dp_test_param['atom_label']

  print ('DP_TEST'.center(80, '*'), flush=True)
  cp2k_frc_file = dp_test_param['cp2k_frc_file']
  cp2k_pos_file = dp_test_param['cp2k_pos_file']
  cp2k_cell_file = dp_test_param['cp2k_cell_file']
  dpff_file = dp_test_param['dpff_file']
  lmp_exe = dp_test_param['lmp_exe']
  lmp_mpi_num = dp_test_param['lmp_mpi_num']
  lmp_omp_num = dp_test_param['lmp_omp_num']

  print ('Do test for supervised learning type', flush=True)

  supervised_test(cp2k_pos_file, cp2k_cell_file, cp2k_frc_file, lmp_exe, lmp_mpi_num, lmp_omp_num, dpff_file, atom_label, work_dir)
