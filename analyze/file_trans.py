#! /usr/env/bin python

import linecache
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.analyze import check_analyze
from CP2K_kit.deepff import gen_lammps_task

def coord2lmp(transd_file, box_file, atom_label, work_dir, file_name):

  '''
  coord2lmp: transform xyz file to lammps data file

  Args:
    transd_file: string
      transd_file is the name of transformed pdb file.
    box_file: string
      box_file is the name of box file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    lmp_file_name: string
      lmp_file_name is the name of transformed lammps data file.
  '''

  tri_cell_vec, atoms, x, y, z = gen_lammps_task.get_box_coord(box_file, transd_file)
  atoms_type_index = []
  for i in atoms:
    atoms_type_index.append(int(data_op.get_dic_keys(atom_label, i)[0]))
  gen_lammps_task.gen_data_file(tri_cell_vec, atoms_type_index, x, y, z, work_dir, file_name)

  lmp_file_name = ''.join((work_dir, '/', file_name))

  return lmp_file_name

def xyz2pdb(transd_file, work_dir, file_name):

  '''
  xyz2pdb: transform pdb file to xyz file

  Args:
    transd_file: string
      transd_file is the name of transformed pdb file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    pdb_file_name: string
      pdb_file_name is the name of transformed pdb file.
  '''

  line = linecache.getline(transd_file, 1)
  atoms_num = int(line.strip('\n'))

  pdb_file_name = ''.join((work_dir, '/', file_name))
  pdb_file = open(pdb_file_name, 'w')

  for i in range(atoms_num):
    line = linecache.getline(transd_file, i+2+1)
    line_split = data_op.split_str(line, ' ', '\n')
    atom = line_split[0]
    atom_label = atom+'R'
    x = float(line_split[1])
    y = float(line_split[2])
    z = float(line_split[3])
    pdb_file.write('%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n' \
                   %('ATOM', i+1, atom_label, 'RES', 'A', 1, x, y, z, 1.00, 0.00, atom))

  linecache.clearcache()
  pdb_file.close()

  return pdb_file_name

def pdb2xyz(transd_file, pre_base, end_base, block_pre_base, block_end_base, time_step, print_freq, work_dir, file_name):

  '''
  pdb2xyz: transform xyz file to pdb file

  Args:
    transd_file: string
      transd_file is the name of transformed xyz file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
    file_name: string
      file_name is the name of generated file.
  Returns:
    xyz_file_name: string
      xyz_file_name is the name of transformed xyz file.
  '''

  line_num = len(open(transd_file).readlines())
  atoms_num = 0
  for i in range(line_num-pre_base-block_pre_base):
    line = linecache.getline(transd_file, i+1+pre_base+block_pre_base)
    line_split = data_op.split_str(line, ' ')
    if ( line_split[0] == 'ATOM' ):
      atoms_num = atoms_num+1
    else:
      break

  frames_num = int((line_num-pre_base-end_base)/(block_pre_base+block_end_base+atoms_num))
  xyz_file_name = ''.join((work_dir, '/', file_name))
  xyz_file = open(xyz_file_name, 'w')

  for i in range(frames_num):
    xyz_file.write('%d\n'%(atoms_num))
    xyz_file.write('%s%9d%s%13.3f%s%21.10f\n' %(' i =', i*print_freq, ', time =', i*time_step*print_freq, ', E =', 0.0))
    for j in range(atoms_num):
      line = linecache.getline(transd_file, i*(block_pre_base+block_end_base+atoms_num)+j+1+pre_base+block_pre_base)
      line_split = data_op.split_str(line, ' ', '\n')
      if ( len(line_split) == 11 ):
        xyz_file.write('%3s%21.10s%20.10s%20.10s\n' %(line_split[10].strip('\n'), line_split[5], line_split[6], line_split[7]))
      elif ( len(line_split) == 12 ):
        xyz_file.write('%3s%21.10s%20.10s%20.10s\n' %(line_split[11].strip('\n'), line_split[6], line_split[7], line_split[8]))

  linecache.clearcache()
  xyz_file.close()

  return xyz_file_name

def file_trans_run(file_trans_param, work_dir):

  '''
  file_trans_run: kernel function to transform file with different format.

  Args:
    file_trans_param: dictionary
      file_trans_param contains keywords used in transforming file.
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    none
  '''

  file_trans_param = check_analyze.check_file_trans_inp(file_trans_param)

  transd_file = file_trans_param['transd_file']
  trans_type = file_trans_param['trans_type']

  print ('FILE_TRANS'.center(80, '*'), flush=True)

  if ( trans_type == 'pdb2xyz' ):
    str_print = 'The %s pdb file will be transfered to xyz type file' %(transd_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    pre_base = file_trans_param['pre_base']
    block_pre_base = file_trans_param['block_pre_base']
    block_end_base = file_trans_param['block_end_base']
    end_base = file_trans_param['end_base']
    time_step = file_trans_param['time_step']
    print_freq = file_trans_param['print_freq']
    xyz_file_name = pdb2xyz(transd_file, pre_base, end_base, block_pre_base, \
                            block_end_base, time_step, print_freq, work_dir, 'coord.xyz')

    str_print = 'The xyz type file is written in %s' %(xyz_file_name)
    print (data_op.str_wrap(str_print, 80), flush=True)

  elif ( trans_type == 'xyz2pdb' ):
    str_print = 'The %s xyz file will be transfered to pdb type file' %(transd_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    print ('The transfered pdb file is a crude file, please revise it in detail!', flush=True)

    pdb_file_name  = xyz2pdb(transd_file, work_dir, 'coord.pdb')

    str_print = 'The pdb type file is written in %s' %(pdb_file_name)
    print (data_op.str_wrap(str_print, 80), flush=True)

  elif ( trans_type == 'coord2lmp' ):
    str_print = 'The %s coord file will be transfered to lammps data type file' %(transd_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    box_file = file_trans_param['box_file']
    atom_label = file_trans_param['atom_label']
    lmp_file_name  = coord2lmp(transd_file, box_file, atom_label, work_dir, 'data.lmp')

    str_print = 'The lammps data type file is written in %s' %(lmp_file_name)
    print (data_op.str_wrap(str_print, 80), flush=True)
