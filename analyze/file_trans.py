#! /usr/env/bin python

import linecache
from CP2K_kit.tools import log_info
from CP2K_kit.tools import list_dic_op

def xyz2pdb(file_name, work_dir):

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

  print ('The transfered pdb file is a crude file, please revise it in detail!', flush=True)

  line = linecache.getline(file_name, 1)
  atoms_num = int(line.strip('\n'))

  pdb_file_name = ''.join((work_dir, '/file.pdb'))
  pdb_file = open(pdb_file_name, 'w')

  for i in range(atoms_num):
    line = linecache.getline(file_name, i+2+1)
    line_split = list_dic_op.str_split(line, ' ')
    atom = line_split[0]
    atom_label = atom+'R'
    x = float(line_split[1])
    y = float(line_split[2])
    z = float(line_split[3].strip('\n'))
    pdb_file.write('%-6s%5d  %-3s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n' \
                   %('ATOM', i+1, atom_label, 'RES', 'A', 1, x, y, z, 1.00, 0.00, atom))

  pdb_file.close()

def pdb2xyz(file_name, work_dir):

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

  line_num = len(open(file_name).readlines())
  for i in range(line_num):
    line = linecache.getline(file_name, i+1)
    line_split = list_dic_op.str_split(line, ' ')
    if ( line_split[0] == 'ATOM' ):
      line_start = i
      break

  xyz_file_name = ''.join((work_dir, '/file.xyz'))
  xyz_file = open(xyz_file_name, 'w')

  total_atom_num = 0
  for i in range(line_num-line_start):
    line = linecache.getline(file_name, i+line_start+1)
    line_split = list_dic_op.str_split(line.strip('\n'), ' ')
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

  if ( 'file_name' in file_trans_param.keys() ):
    file_name = file_trans_param['file_name']
  else:
    log_info.log_error('No file need to be transfered, please set analzye/file_trans/file_name')
    exit()

  if ( 'trans_type' in file_trans_param.keys() ):
    trans_type = file_trans_param['trans_type']
  else:
    log_info.log_error('No transfer type found, please set analyze/file_trans/trans_type')
    exit()

  if ( trans_type == 'pdb2xyz' ):
    pdb2xyz(file_name, work_dir)
  elif ( trans_type == 'xyz2pdb' ):
    xyz2pdb(file_name, work_dir)
