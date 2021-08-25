#! /usr/env/bin python

import linecache
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.analyze import check_analyze

def xyz2pdb(transd_file, work_dir, file_name):

  '''
  xyz2pdb : transfer pdb file to xyz file

  Args :
    file_name : string
      file_name is the name of transformed pdb file.
    work_dir : string
      work_dir is the working directory
  Returns :
    none
  '''

  line = linecache.getline(transd_file, 1)
  atoms_num = int(line.strip('\n'))

  pdb_file_name = ''.join((work_dir, '/', file_name))
  pdb_file = open(pdb_file_name, 'w')

  for i in range(atoms_num):
    line = linecache.getline(transd_file, i+2+1)
    line_split = data_op.str_split(line, ' ')
    atom = line_split[0]
    atom_label = atom+'R'
    x = float(line_split[1])
    y = float(line_split[2])
    z = float(line_split[3].strip('\n'))
    pdb_file.write('%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n' \
                   %('ATOM', i+1, atom_label, 'RES', 'A', 1, x, y, z, 1.00, 0.00, atom))

  pdb_file.close()

  return pdb_file_name

def pdb2xyz(transd_file, work_dir, file_name):

  '''
  pdb2xyz : transfer xyz file to pdb file

  Args :
    file_name : string
      file_name is the name of transformed xyz file.
    work_dir : string
      work_dir is the working directory
  Returns :
    none
  '''

  line_num = len(open(transd_file).readlines())
  for i in range(line_num):
    line = linecache.getline(transd_file, i+1)
    line_split = data_op.str_split(line, ' ')
    if ( line_split[0] == 'ATOM' ):
      line_start = i
      break

  xyz_file_name = ''.join((work_dir, '/', file_name))
  xyz_file = open(xyz_file_name, 'w')

  total_atom_num = 0
  for i in range(line_num-line_start):
    line = linecache.getline(transd_file, i+line_start+1)
    line_split = data_op.str_split(line.strip('\n'), ' ')
    if ( len(line_split) == 11 ):
      total_atom_num = total_atom_num + 1
      xyz_file.write('%-3s%8s%8s%8s\n' %(line_split[10].strip('\n'), line_split[5], line_split[6], line_split[7]))
    elif ( len(line_split) == 12 ):
      total_atom_num = total_atom_num + 1
      xyz_file.write('%-3s%8s%8s%8s\n' %(line_split[11].strip('\n'), line_split[6], line_split[7], line_split[8]))

  xyz_file.close()

  xyz_file = open(xyz_file_name, 'r+')
  content = xyz_file.read()
  xyz_file.seek(0, 0)
  xyz_file.write('%d\n'%(total_atom_num) + '\n' + content)
  xyz_file.close()

  return xyz_file_name

def file_trans_run(file_trans_param, work_dir):

  '''
  file_trans_run : kernel function to transform file with different format

  Args :
    file_trans_param : dictionary
      file_trans_param contains keywords used in transforming file.
    work_dir : string
      work_dir is working directory of CP2K_kit.
  Returns :
    none
  '''

  file_trans_param = check_analyze.check_file_trans_inp(file_trans_param)

  transd_file = file_trans_param['transd_file']
  trans_type = file_trans_param['trans_type']

  print ('FILE_TRANS'.center(80, '*'), flush=True)

  if ( trans_type == 'pdb2xyz' ):
    str_print = 'The %s pdb file will be transfered to xyz type file' %(transd_file)
    print (data_op.str_wrap(str_print, 80), flush=True)

    xyz_file_name = pdb2xyz(transd_file, work_dir, 'coord.xyz')

    str_print = 'The xyz type file is written in %s' %(xyz_file_name)
    print (data_op.str_wrap(str_print, 80), flush=True)

  elif ( trans_type == 'xyz2pdb' ):
    str_print = 'The %s xyz file will be transfered to pdb type file' %(transd_file)
    print (data_op.str_wrap(str_print, 80), flush=True)
    print ('The transfered pdb file is a crude file, please revise it in detail!', flush=True)

    pdb_file_name  = xyz2pdb(transd_file, work_dir, 'coord.pdb')

    str_print = 'The pdb type file is written in %s' %(pdb_file_name)
    print (data_op.str_wrap(str_print, 80), flush=True)
