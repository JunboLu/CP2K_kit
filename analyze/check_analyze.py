#! /use/env/bin python

import os
import copy
from collections import OrderedDict
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.tools import traj_info

def check_center_inp(center_dic):

  '''
  check_center_inp : check the input file of center.

  Args :
    center_dic : dictionary
      center_dic contains the parameter for center.
  Returns :
    new_center_dic : dictionary
      new_center_dic is the revised center_dic
  '''

  #As we use pop, so we copy the dic.
  new_center_dic = copy.deepcopy(center_dic)

  if ( 'center_type' in new_center_dic.keys() ):
    center_type = new_center_dic['center_type']
    if ( center_type == 'center_box' or center_type == 'center_image' ):
      pass
    else:
      log_info.log_error('Input error: only center_box and center_image are supported, please check or set analyze/center/type ')
      exit()
  else:
    log_info.log_error('Input error: no center type, please set analyze/center/type')
    exit()

  if ( new_center_dic['center_type'] == 'center_image' ):
    if ( 'center_atom_id' in new_center_dic.keys() ):
      center_id = new_center_dic['center_atom_id']
      if ( data_op.eval_str(center_id) == 1 ):
        new_center_dic['center_atom_id'] = int(center_id)
      else:
        log_info.log_error('Input error: center atom id should be integer, please check or set analyze/center/center_id')
        exit()
    else:
      log_info.log_error('Input error: no center atom id for center_image, please set analyze/center/center_id')
      exit()

  if ( 'traj_coord_file' in new_center_dic.keys() ):
    traj_coord_file = new_center_dic['traj_coord_file']
    if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
      new_center_dic['traj_file'] = os.path.abspath(traj_coord_file)
    else:
      log_info.log_error('%s file does not exist' %(traj_coord_file))
      exit()
  else:
    log_info.log_error('Input error: no coordination trajectory file, please set analyze/center/traj_coord_file')
    exit()

  if ( 'box' in new_center_dic.keys() ):
    A_exist = 'A' in new_center_dic['box'].keys()
    B_exist = 'B' in new_center_dic['box'].keys()
    C_exist = 'C' in new_center_dic['box'].keys()
  else:
    log_info.log_error('Input error: no box, please set analyze/center/box')
    exit()

  if ( A_exist and B_exist and C_exist ):
    box_A = new_center_dic['box']['A']
    box_B = new_center_dic['box']['B']
    box_C = new_center_dic['box']['C']
  else:
    log_info.log_error('Input error: box setting error, please check analyze/center/box')
    exit()

  if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
    new_center_dic['box']['A'] = [float(x) for x in box_A]
  else:
    log_info.log_error('Input error: A vector of box wrong, please check analyze/center/box/A')
    exit()

  if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
    new_center_dic['box']['B'] = [float(x) for x in box_B]
  else:
    log_info.log_error('Input error: B vector of box wrong, please check analyze/center/box/B')
    exit()

  if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
    new_center_dic['box']['C'] = [float(x) for x in box_C]
  else:
    log_info.log_error('Input error: C vector of box wrong, please check analyze/center/box/C')
    exit()

  if ( 'connect0' in new_center_dic.keys() ):
    group_atom = []
    atom_id = []
    group_num = 0
    for i in new_center_dic['connect'].keys():
      if ( 'group' in i ):
        group_num = group_num+1
        if ( 'atom_id' in new_center_dic['connect'][i].keys() ):
          atom_id_i = data_op.get_id_list(new_center_dic['connect'][i]['atom_id'])
          atom_id.append(atom_id_i)
        else:
          log_info.log_error('Input error: no atom id, please set analyze/center/connect/group/atom_id')
          exit()
        if ( 'group_atom' in new_center_dic['connect'][i].keys() ):
          group_atom_i = new_center_dic['connect'][i]['group_atom']
          if ( isinstance(group_atom_i, list)):
            if ( all(data_op.eval_str(x) == 0 for x in group_atom_i) ):
              group_atom.append(group_atom_i)
            else:
              log_info.log_error('Input error: group atoms wrong, please check or reset analyze/center/connect/group/group_atom')
              exit()
          else:
            group_atom.append([group_atom_i])
        else:
          log_info.log_error('Input error: no group atoms, please set analyze/center/connect/group/group_atom')
          exit()

    for i in center_dic['connect0'].keys():
      new_center_dic['connect0'].pop(i)

  new_center_dic['connect0']['atom_id'] = atom_id
  new_center_dic['connect0']['group_atom'] = group_atom

  return new_center_dic

def check_diffusion_inp(diffusion_dic):

  '''
  check_diffusion_inp : check the input of diffusion.

  Args :
    diffusion_dic : dictionary
      diffusion_dic contains parameters for diffusion
  Returns :
    diffusion_dic: dictionary
      diffusion_dic is the revised diffusion_dic
  '''

  #new_diffusion_dic = copy.deepcopy(diffusion_dic)

  if ( 'method' in diffusion_dic.keys() ):
    method = diffusion_dic['method']
    if ( method == 'einstein_sum' or method == 'green_kubo' ):
      pass
    else:
      log_info.log_error('Input error: only einstein_sum or green_kubo are supported for diffusion calculation')
      exit()
  else:
    diffusion_dic['method'] = 'einstein_sum'

  method = diffusion_dic['method']
  if ( method == 'einstein_sum' ):
    if ( 'traj_coord_file' in diffusion_dic.keys() ):
      traj_coord_file = diffusion_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        diffusion_dic['traj_coord_file'] = os.path.abspath(traj_coord_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_coord_file, 'coord')
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
        exit()
    else:
      log_info.log_error('Input error: no traj_file, please set analyze/diffusion/traj_coord_file')
      exit()

    if ( 'remove_com' in diffusion_dic.keys() ):
      remove_com = data_op.str_to_bool(diffusion_dic['remove_com'])
      if ( isinstance(remove_com, bool) ):
        diffusion_dic['remove_com'] = remove_com
      else:
        log_info.log_error('Input error: remove_com must be bool, please check or reset analyze/diffusion/remove_com')
    else:
      diffusion_dic['remove_com'] = True

  elif ( method == 'green_kubo' ):
    if ( 'traj_vel_file' in diffusion_dic.keys() ):
      traj_vel_file = diffusion_dic['traj_vel_file']
      if ( os.path.exists(os.path.abspath(traj_vel_file)) ):
        diffusion_dic['traj_file'] = os.path.abspath(traj_vel_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_vel_file, 'vel')
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_vel_file))
        exit()
    else:
      log_info.log_error('Input error: no traj_file, please set analyze/diffusion/traj_vel_file')
      exit()

  if ( 'atom_id' in diffusion_dic.keys() ):
    atom_id = data_op.get_id_list(diffusion_dic['atom_id'])
    diffusion_dic['atom_id'] = atom_id
  else:
    log_info.log_error('Input error: no atom_id, please set analyze/diffusion/atom_id')
    exit()

  if ( 'init_step' in diffusion_dic.keys() ):
    init_step = diffusion_dic['init_step']
    if ( data_op.eval_str(init_step) == 1 ):
      diffusion_dic['init_step'] = int(init_step)
    else:
      log_info.log_error('Input error: init_step wrong, please check or set analyze/diffusion/init_step')
      exit()
  else:
    diffusion_dic['init_step'] = start_frame_id

  if ( 'end_step' in diffusion_dic.keys() ):
    end_step = diffusion_dic['end_step']
    if ( data_op.eval_str(end_step) == 1 ):
      diffusion_dic['end_step'] = int(end_step)
    else:
      log_info.log_error('Input error: end_step wrong, please check or set analyze/diffusion/end_step')
      exit()
  else:
    diffusion_dic['end_step'] = end_frame_id

  if ( 'max_frame_corr' in diffusion_dic.keys() ):
    max_frame_corr = diffusion_dic['max_frame_corr']
    if ( data_op.eval_str(max_frame_corr) == 1 ):
      if ( int(max_frame_corr) > int(frames_num/2) ):
        log_info.log_error('Input error: max_frame_corr should be less than frames_num/2, please check or reset analyze/diffusion/max_frame_corr')
        exit()
      else:
        diffusion_dic['max_frame_corr'] = int(max_frame_corr)
    else:
      log_info.log_error('Input error: max_frame_corr should be integer, please check or set analyze/diffusion/max_frame_corr')
      exit()
  else:
    diffusion_dic['max_frame_corr'] = int(frames_num/2)

  return diffusion_dic

def check_file_trans_inp(file_trans_dic):

  '''
  check_file_trans_inp : check the input of file_trans.

  Args :
    file_trans_dic : dictionary
      file_trans_dic contains parameters for file_trans.
  Returns :
    file_trans_dic: dictionary
      file_trans_dic is the revised file_trans_dic
  '''

  if ( 'transd_file' in file_trans_dic.keys() ):
    transd_file = file_trans_dic['transd_file']
    if ( os.path.exists(os.path.abspath(transd_file)) ):
      file_trans_dic['transd_file'] = os.path.abspath(transd_file)
    else:
      log_info.log_error('Input error: %s does not exist' %(transd_file))
      exit()
  else:
    log_info.log_error('Input error: no transfered file, please set analzye/file_trans/transd_file')
    exit()

  if ( 'trans_type' in file_trans_dic.keys() ):
    trans_type = file_trans_dic['trans_type']
    if ( trans_type == 'pdb2xyz' or trans_type == 'xyz2pdb' ):
      pass
    else:
      log_info.log_error('Input error: only pbd2xyz and xyz2pdb are supported, please check or reset analyze/file_trans/trans_type')
      exit()
  else:
    log_info.log_error('Input error: no transfer type, please set analyze/file_trans/trans_type')
    exit()

  return file_trans_dic

def check_geometry_inp(geometry_dic):

  '''
  check_geometry_inp : check the input of geometry.

  Args :
    geometry_dic : dictionary
      geometry_dic contains parameters for geometry.
  Returns :
    geometry_dic: dictionary
      geometry_dic is the revised geometry_dic
  '''

  if ( 'coord_num' in geometry_dic ):
    coord_num_dic = geometry_dic['coord_num']

    if ( 'traj_coord_file' in coord_num_dic.keys() ):
      traj_coord_file = coord_num_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        geometry_dic['coord_num']['traj_coord_file'] = os.path.abspath(traj_coord_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(os.path.abspath(traj_coord_file), 'coord')
      else:
        log_info.log_error('Input error: %s does not exist' %(traj_coord_file))
        exit()
    else:
      log_info.log_error('Input error: no coordination trajectory file, please set analyze/geometry/coord_num/traj_coord_file')
      exit()

    if ( 'init_step' in coord_num_dic.keys() ):
      init_step = coord_num_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        geometry_dic['coord_num']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step should be integer, please check or reset analyze/geometry/coord_num/init_step')
        exit()
    else:
      geometry_dic['coord_num']['init_step'] = start_frame_id

    if ( 'end_step' in coord_num_dic.keys() ):
      end_step = coord_num_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        geometry_dic['coord_num']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step should be integer, please check or reset analyze/geometry/coord_num/end_step')
        exit()
    else:
      geometry_dic['coord_num']['end_step'] = end_frame_id

    if ( 'r_cut' in coord_num_dic.keys() ):
      r_cut = coord_num_dic['r_cut']
      if ( data_op.eval_str(r_cut) == 1 or data_op.eval_str(r_cut) ==2 ):
        geometry_dic['coord_num']['r_cut'] = float(r_cut)
      else:
        log_info.log_error('Input error: r_cut must be float, please check or reset analyze/geometry/coord_num/r_cut')
    else:
      geometry_dic['coord_num']['r_cut'] = 6.0

    if ( 'box' in coord_num_dic.keys() ):
      A_exist = 'A' in coord_num_dic['box'].keys()
      B_exist = 'B' in coord_num_dic['box'].keys()
      C_exist = 'C' in coord_num_dic['box'].keys()
    else:
      log_info.log_error('Input error: no box, please set analyze/geometry/coord_num/box')
      exit()

    if ( A_exist and B_exist and C_exist ):
      box_A = coord_num_dic['box']['A']
      box_B = coord_num_dic['box']['B']
      box_C = coord_num_dic['box']['C']
    else:
      log_info.log_error('Input error: box setting error, please check analyze/geometry/coord_num/box')
      exit()

    if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
      geometry_dic['coord_num']['box']['A'] = [float(x) for x in box_A]
    else:
      log_info.log_error('Input error: A vector of box wrong, please check analyze/geometry/coord_num//box/A')
      exit()

    if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
      geometry_dic['coord_num']['box']['B'] = [float(x) for x in box_B]
    else:
      log_info.log_error('Input error: B vector of box wrong, please check analyze/geometry/coord_num/box/B')
      exit()

    if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
      geometry_dic['coord_num']['box']['C'] = [float(x) for x in box_C]
    else:
      log_info.log_error('Input error: C vector of box wrong, please check analyze/geometry/coord_num/box/C')
      exit()

    return geometry_dic

  elif ( 'bond_length' in geometry_dic ):
    bond_length_dic = geometry_dic['bond_length']

    if ( 'traj_coord_file' in bond_length_dic.keys() ):
      traj_coord_file = bond_length_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        geometry_dic['bond_length']['traj_coord_file'] = os.path.abspath(traj_coord_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(os.path.abspath(traj_coord_file), 'coord')
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
        exit()
    else:
      log_info.log_error('Input error: no coordination trajectory file, please set analyze/geometry/bond_length/traj_coord_file')
      exit()

    if ( 'atom_pair' in bond_length_dic.keys() ):
      atom_pair = bond_length_dic['atom_pair']
      if ( len(atom_pair) == 2 and all(data_op.eval_str(x) == 1 for x in atom_pair) ):
        geometry_dic['bond_length']['atom_pair'] = [int (x) for x in atom_pair]
      else:
        log_info.log_error('Input error: atom_pair should be 2 integer, please check or reset analyze/geometry/bond_length/atom_pair')
        exit()
    else:
      log_info.log_error('Input error: no atom_pair, please set analyze/geometry/bond_length/atom_pair')
      exit()

    if ( 'init_step' in bond_length_dic.keys() ):
      init_step = bond_length_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        geometry_dic['bond_length']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step should be integer, please check or reset analyze/geometry/bond_length/init_step')
        exit()
    else:
      geometry_dic['bond_length']['init_step'] = start_frame_id

    if ( 'end_step' in bond_length_dic.keys() ):
      end_step = bond_length_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        geometry_dic['bond_length']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step should be integer, please check or reset analyze/geometry/bond_length/end_step')
        exit()
    else:
      geometry_dic['bond_length']['end_step'] = end_frame_id

    if ( 'box' in bond_length_dic.keys() ):
      A_exist = 'A' in bond_length_dic['box'].keys()
      B_exist = 'B' in bond_length_dic['box'].keys()
      C_exist = 'C' in bond_length_dic['box'].keys()
    else:
      log_info.log_error('Input error: no box, please set analyze/geometry/bond_length/box')
      exit()

    if ( A_exist and B_exist and C_exist ):
      box_A = bond_length_dic['box']['A']
      box_B = bond_length_dic['box']['B']
      box_C = bond_length_dic['box']['C']
    else:
      log_info.log_error('Input error: box setting error, please check analyze/geometry/bond_length/box')
      exit()

    if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
      geometry_dic['bond_length']['box']['A'] = [float(x) for x in box_A]
    else:
      log_info.log_error('Input error: A vector of box wrong, please check analyze/geometry/bond_length//box/A')
      exit()

    if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
      geometry_dic['bond_length']['box']['B'] = [float(x) for x in box_B]
    else:
      log_info.log_error('Input error: B vector of box wrong, please check analyze/geometry/bond_length/box/B')
      exit()

    if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
      geometry_dic['bond_length']['box']['C'] = [float(x) for x in box_C]
    else:
      log_info.log_error('Input error: C vector of box wrong, please check analyze/geometry/bond_length/box/C')
      exit()

    return geometry_dic

  elif ( 'bond_angle' in geometry_dic ):
    bond_angle_dic = geometry_dic['bond_angle']

    if ( 'traj_coord_file' in bond_angle_dic.keys() ):
      traj_coord_file = bond_angle_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        geometry_dic['bond_angle']['traj_coord_file'] = os.path.abspath(traj_coord_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_coord_file, 'coord')
      else:
        log_info.log_error('Input error: %s does not exist' %(traj_coord_file))
    else:
      log_info.log_error('Input error: no coordination trajectory file, please set analyze/geometry/bond_angle/traj_coord_file')
      exit()

    if ( 'atom_pair' in bond_angle_dic.keys() ):
      atom_pair = bond_angle_dic['atom_pair']
      if ( len(atom_pair) == 3 and all(data_op.eval_str(x) == 1 for x in atom_pair) ):
        geometry_dic['bond_angle']['atom_pair'] = [int(x) for x in atom_pair]
      else:
        log_info.log_error('Input error: atom_pair should be 3 integers, please check or reset analyze/geometry/bond_angle/atom_pair')
        exit()
    else:
      log_info.log_error('Input error: no atom_pair, please set analyze/geometry/bond_angle/atom_pair')
      exit()

    if ( 'init_step' in bond_angle_dic.keys() ):
      init_step = bond_angle_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        geometry_dic['bond_angle']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step shoule be integer, please check or reset analyze/geometry/bond_angle/init_step')
        exit()
    else:
      geometry_dic['bond_angle']['init_step'] = start_frame_id

    if ( 'end_step' in bond_angle_dic.keys() ):
      end_step = bond_angle_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        geometry_dic['bond_angle']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step shoule be integer, please check or reset analyze/geometry/bond_angle/end_step')
        exit()
    else:
      geometry_dic['bond_angle']['end_step'] = end_frame_id

    return geometry_dic

  elif ( 'first_shell' in geometry_dic ):
    first_shell_dic = geometry_dic['first_shell']

    if ( 'traj_coord_file' in first_shell_dic.keys() ):
      traj_coord_file = first_shell_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        geometry_dic['first_shell']['traj_coord_file'] = os.path.abspath(traj_coord_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(os.path.abspath(traj_coord_file), 'coord')
      else:
        log_info.log_error('Input error: %s does not exist' %(traj_coord_file))
        exit()

    if ( 'atom_type_pair' in first_shell_dic.keys() ):
      atom_type_pair = first_shell_dic['atom_type_pair']
      if ( len(atom_type_pair) == 2 and all(data_op.eval_str(x) == 0 for x in atom_type_pair) ):
        pass
      else:
        log_info.log_error('Input error: atom_type_pair should be 2 string, please check or reset analyze/geometry/first_shell/atom_type_pair')
        exit()
    else:
      log_info.log_error('Input error: no atom_type_pair, please set analyze/geometry/first_shell/atom_type_pair')
      exit()

    if ( 'first_shell_dist' in first_shell_dic.keys() ):
      first_shell_dist = first_shell_dic['first_shell_dist']
      if ( data_op.eval_str(first_shell_dist) == 2 ):
        geometry_dic['first_shell']['first_shell_dist'] = float(first_shell_dist)
      else:
        log_info.log_error('Input error: first_shell_dist should be float, please check or reset analyze/geometry/first_shell/first_shell_dist')
        exit()
    else:
      log_info.log_error('Input error: no first_shell_dist, please set analyze/geometry/first_shell/first_shell_dist')
      exit()

    if ( 'dist_conv' in first_shell_dic.keys() ):
      dist_conv = first_shell_dic['dist_conv']
      if ( data_op.eval_str(dist_conv) == 2 ):
        geometry_dic['first_shell']['dist_conv'] = float(dist_conv)
      else:
        log_info.log_error('Input error: dist_conv should be float, please check or reset analyze/geometry/first_shell/dist_conv')
        exit()
    else:
      geometry_dic['first_shell']['dist_conv'] = 0.1

    if ( 'init_step' in first_shell_dic.keys() ):
      init_step = first_shell_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        geometry_dic['first_shell']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step should be integer, please check or reset analyze/geometry/first_shell/init_step')
    else:
      geometry_dic['first_shell']['init_step'] = start_frame_id

    if ( 'end_step' in first_shell_dic.keys() ):
      end_step = first_shell_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        geometry_dic['first_shell']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step should be integer, please check or reset analyze/geometry/first_shell/end_step')
    else:
      geometry_dic['first_shell']['end_step'] = end_frame_id

    if ( 'box' in first_shell_dic.keys() ):
      A_exist = 'A' in first_shell_dic['box'].keys()
      B_exist = 'B' in first_shell_dic['box'].keys()
      C_exist = 'C' in first_shell_dic['box'].keys()
    else:
      log_info.log_error('Input error: no box, please set analyze/geometry/first_shell/box')
      exit()

    if ( A_exist and B_exist and C_exist ):
      box_A = first_shell_dic['box']['A']
      box_B = first_shell_dic['box']['B']
      box_C = first_shell_dic['box']['C']
    else:
      log_info.log_error('Input error: box setting error, please check analyze/geometry/first_shell/box')
      exit()

    if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
      geometry_dic['first_shell']['box']['A'] = [float(x) for x in box_A]
    else:
      log_info.log_error('Input error: A vector of box wrong, please check analyze/geometry/first_shell//box/A')
      exit()

    if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
      geometry_dic['first_shell']['box']['B'] = [float(x) for x in box_B]
    else:
      log_info.log_error('Input error: B vector of box wrong, please check analyze/geometry/first_shell/box/B')
      exit()

    if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
      geometry_dic['first_shell']['box']['C'] = [float(x) for x in box_C]
    else:
      log_info.log_error('Input error: C vector of box wrong, please check analyze/geometry/first_shell/box/C')
      exit()

    return geometry_dic

  elif ( 'choose_structure' in geometry_dic ):
    choose_str_dic = geometry_dic['choose_structure']

    if ( 'file_type' in choose_str_dic.keys() ):
      file_type = choose_str_dic['file_type']
      valid_file_type = ['coord', 'vel', 'frc']
      if ( file_type in valid_file_type ):
        pass
      else:
        log_info.log_error('Input error: only coord, vel, and frc are supported for file_type,\
                            please check or reset analyze/geometry/choose_structure/file_type')
        exit()
    else:
      log_info.log_error('Input error: no file type, please set analyze/geometry/choose_structure/file_type')
      exit()
    file_type = choose_str_dic['file_type']

    if ( 'traj_file' in choose_str_dic.keys() ):
      traj_file = choose_str_dic['traj_file']
      if ( os.path.exists(os.path.abspath(traj_file)) ):
        geometry_dic['choose_structure']['traj_file'] = os.path.abspath(traj_file)
        atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
        traj_info.get_traj_info(traj_file, 'coord')
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_file))
        exit()
    else:
      log_info.log_error('Input error: no coordination trajectory file, please set analyze/geometry/choose_structure/traj_file')
      exit()

    if ( 'init_step' in choose_str_dic.keys() ):
      init_step = choose_str_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        geometry_dic['choose_structure']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step should be integer, please check or reset analyze/geometry/choose_structure/init_step')
    else:
      geometry_dic['choose_structure']['init_step'] = start_frame_id

    if ( 'end_step' in choose_str_dic.keys() ):
      end_step = choose_str_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        geometry_dic['choose_structure']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step should be integer, please check or reset analyze/geometry/choose_structure/end_step')
    else:
      geometry_dic['choose_structure']['end_step'] = end_frame_id

    if ( 'atom_id' in choose_str_dic.keys() ):
      geometry_dic['choose_structure']['atom_id'] = data_op.get_id_list(choose_str_dic['atom_id'])
    else:
      log_info.log_error('Input error: no atom_id, please set analyze/geometry/choose_structure/atom_id')
      exit()

    return geometry_dic

  elif ( 'order_structure' in geometry_dic.keys() ):
    new_geometry_dic = copy.deepcopy(geometry_dic)
    order_str_dic = new_geometry_dic['order_structure']

    if ( 'traj_coord_file' in order_str_dic.keys() ):
      traj_coord_file = order_str_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        new_geometry_dic['order_structure']['traj_coord_file'] = os.path.abspath(traj_coord_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
        exit()
    else:
      log_info.log_error('Input error: no coordination trajectory file, please set analyze/geometry/order_structure/traj_coord_file')
      exit()

    if ( 'box' in order_str_dic.keys() ):
      A_exist = 'A' in order_str_dic['box'].keys()
      B_exist = 'B' in order_str_dic['box'].keys()
      C_exist = 'C' in order_str_dic['box'].keys()
    else:
      log_info.log_error('Input error: no box, please set analyze/geometry/order_structure/box')
      exit()

    if ( A_exist and B_exist and C_exist ):
      box_A = order_str_dic['box']['A']
      box_B = order_str_dic['box']['B']
      box_C = order_str_dic['box']['C']
    else:
      log_info.log_error('Input error: box setting error, please check analyze/geometry/order_structure/box')
      exit()

    if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
      new_geometry_dic['order_structure']['box']['A'] = [float(x) for x in box_A]
    else:
      log_info.log_error('Input error: A vector of box wrong, please check analyze/geometry/order_structure//box/A')
      exit()

    if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
      new_geometry_dic['order_structure']['box']['B'] = [float(x) for x in box_B]
    else:
      log_info.log_error('Input error: B vector of box wrong, please check analyze/geometry/order_structure/box/B')
      exit()

    if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
      new_geometry_dic['order_structure']['box']['C'] = [float(x) for x in box_C]
    else:
      log_info.log_error('Input error: C vector of box wrong, please check analyze/geometry/order_structure/box/C')
      exit()

    group_atom = []
    atom_id = []
    group_num = 0
    for i in order_str_dic.keys():
      if ( 'connect' in i ):
        group_num = group_num+1
        if ( 'atom_id' in order_str_dic[i].keys() ):
          atom_id_i = data_op.get_id_list(order_str_dic[i]['atom_id'])
          atom_id.append(atom_id_i)
        else:
          log_info.log_error('Input error: no atom_id, please set analyze/geometry/order_structure/connect/atom_id')
          exit()
        if ( 'group_atom' in order_str_dic[i].keys() ):
          group_atom_i = order_str_dic[i]['group_atom']
          if ( isinstance(group_atom_i, list)):
            if ( all(data_op.eval_str(x) == 0 for x in group_atom_i) ):
              group_atom.append(group_atom_i)
            else:
              log_info.log_error('Input error: group_atom wrong, please check or reset analyze/geometry/order_structure/connect/group_atom')
              exit()
          else:
            group_atom.append([group_atom_i])
        else:
          log_info.log_error('Input error: no group_atom, please set analyze/geometry/order_structure/connect/group_atom')
          exit()

    if ( group_num == 0 ):
      log_info.log_error('Input error: no connect information, please check or reset analyze/geometry/order_structure/connect')
      exit()
    else:
      for i in geometry_dic['order_structure'].keys():
        if ( 'connect' in i ):
          new_geometry_dic['order_structure'].pop(i)

      new_geometry_dic['order_structure']['atom_id'] = atom_id
      new_geometry_dic['order_structure']['group_atom'] = group_atom

    return new_geometry_dic

def check_lmp2cp2k_inp(lmp2cp2k_dic):

  '''
  check_lmp2cp2k_inp : check the input of lmp2cp2k.

  Args :
    lmp2cp2k_dic : dictionary
      lmp2cp2k_dic contains parameters for lmp2cp2k.
  Returns :
    lmp2cp2k_dic: dictionary
      lmp2cp2k_dic is the revised lmp2cp2k_dic.
  '''

  lmp2cp2k_dic = copy.deepcopy(lmp2cp2k_dic)

  if ( 'lmp_log_file' in lmp2cp2k_dic.keys() ):
    lmp_log_file = lmp2cp2k_dic['lmp_log_file']
    if ( os.path.exists(os.path.abspath(lmp_log_file)) ):
      lmp2cp2k_dic['lmp_log_file'] = os.path.abspath(lmp_log_file)
    else:
      log_info.log_error('Input error: %s file does not exist' %(lmp_log_file))
      exit()
  else:
    log_info.log_error('Input error: no lmp_log_file, please set analyze/lmp2cp2k/lmp_log_file')
    exit()

  if ( 'lmp_traj_file' in lmp2cp2k_dic.keys() ):
    lmp_traj_file = lmp2cp2k_dic['lmp_traj_file']
    if ( os.path.exists(os.path.abspath(lmp_traj_file)) ):
      lmp2cp2k_dic['lmp_traj_file'] = os.path.abspath(lmp_traj_file)
    else:
      log_info.log_error('Input error: %s file does not exist' %(lmp_traj_file))
      exit()
  else:
    log_info.log_error('Input error: no lmp_traj_file, please set analyze/lmp2cp2k/lmp_traj_file')
    exit()

  if ( 'atom_label' in lmp2cp2k_dic.keys() ):
    atom_label = lmp2cp2k_dic['atom_label']
    atom_label_dic = OrderedDict()
    for i in range (len(atom_label)):
      lable_split = data_op.str_split(atom_label[i], ':')
      atom_label_dic[int(lable_split[0])] = lable_split[1]
    lmp2cp2k_dic['atom_label'] = atom_label_dic
  else:
    log_info.log_error('Input error: no atom_label, please set analyze/lmp2cp2k/atom_label')
    exit()

  if ( 'time_step' in lmp2cp2k_dic.keys() ):
    time_step = lmp2cp2k_dic['time_step']
    if ( data_op.eval_str(time_step) == 1 or data_op.eval_str(time_step) == 2 ):
      lmp2cp2k_dic['time_step'] = float(time_step)
    else:
      log_info.log_error('Input error: time_step should float, please set analyze/lmp2cp2k/time_step')
      exit()
  else:
    log_info.log_error('Input error: no time_step, please set analyze/lmp2cp2k/time_step')
    exit()

  valid_unit = ['metal', 'real']
  if ( 'lmp_unit' in lmp2cp2k_dic.keys() ):
    lmp_unit = lmp2cp2k_dic['lmp_unit']
    if ( lmp_unit in valid_unit ):
      pass
    else:
      log_info.log_error('Input error: %s is not supported, please check or reset analyze/lmp2cp2k/lmp_unit')
      exit()
  else:
    log_info.log_error('Input error: no lmp_unit, please set analyze/lmp2cp2k/lmp_unit')
    exit()

  if ( 'unwrap' in lmp2cp2k_dic.keys() ):
    unwrap = data_op.str_to_bool(lmp2cp2k_dic['unwrap'])
    if ( isinstance(unwrap, bool) ):
      lmp2cp2k_dic['unwrap'] = unwrap
    else:
      log_info.log_error('Input error: unwrap should be bool, please set analyze/lmp2cp2k/unwrap')
      exit()
  else:
    lmp2cp2k_dic['unwrap'] = False

  if ( 'box' in lmp2cp2k_dic.keys() ):
    A_exist = 'A' in lmp2cp2k_dic['box'].keys()
    B_exist = 'B' in lmp2cp2k_dic['box'].keys()
    C_exist = 'C' in lmp2cp2k_dic['box'].keys()
  else:
    log_info.log_error('Input error: no box, please set analyze/lmp2cp2k/box')
    exit()

  if ( A_exist and B_exist and C_exist ):
    box_A = lmp2cp2k_dic['box']['A']
    box_B = lmp2cp2k_dic['box']['B']
    box_C = lmp2cp2k_dic['box']['C']
  else:
    log_info.log_error('Input error: box setting error, please check analyze/lmp2cp2k/box')
    exit()

  if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
    lmp2cp2k_dic['box']['A'] = [float(x) for x in box_A]
  else:
    log_info.log_error('Input error: A vector of box wrong, please check analyze/lmp2cp2k/box/A')
    exit()

  if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
    lmp2cp2k_dic['box']['B'] = [float(x) for x in box_B]
  else:
    log_info.log_error('Input error: B vector of box wrong, please check analyze/lmp2cp2k/box/B')
    exit()

  if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
    lmp2cp2k_dic['box']['C'] = [float(x) for x in box_C]
  else:
    log_info.log_error('Input error: C vector of box wrong, please check analyze/lmp2cp2k/box/C')
    exit()

  return lmp2cp2k_dic

def check_rdf_inp(rdf_dic):

  '''
  check_lmp2cp2k_inp : check the input of rdf.

  Args :
    rdf_dic : dictionary
      rdf_dic contains parameters for rdf.
  Returns :
    rdf_dic: dictionary
      rdf_dic is the revised rdf_dic.
  '''

  rdf_dic = copy.deepcopy(rdf_dic)

  if ( 'traj_coord_file' in rdf_dic.keys() ):
    traj_coord_file = rdf_dic['traj_coord_file']
    if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
      rdf_dic['traj_coord_file'] = os.path.abspath(traj_coord_file)
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(os.path.abspath(traj_coord_file), 'coord')
    else:
      log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
      exit()
  else:
    log_info.log_error('Input error: no coordination trajectroy file, please set analyze/rdf/traj_coord_file')
    exit()

  if ( 'atom_type_pair' in rdf_dic.keys() ):
    atom_type_pair = rdf_dic['atom_type_pair']
    if ( len(atom_type_pair) == 2 and all(data_op.eval_str(x) == 0 for x in atom_type_pair) ):
      pass
    else:
      log_info.log_error('Input error: atom_type_pair should be 2 string, please check or reset analyze/rdf/atom_type_pair')
      exit()
  else:
    log_info.log_error('Input error: no atom type, please set analyze/rdf/atom_type_pair')
    exit()

  if ( 'init_step' in rdf_dic.keys() ):
    init_step = rdf_dic['init_step']
    if ( data_op.eval_str(init_step) == 1 ):
      rdf_dic['init_step'] = int(init_step)
    else:
      log_info.log_error('Input error: init_step should be integer, please check or reset analyze/rdf/init_step')
  else:
    rdf_dic['init_step'] = start_frame_id

  if ( 'end_step' in rdf_dic.keys() ):
    end_step = rdf_dic['end_step']
    if ( data_op.eval_str(end_step) == 1 ):
      rdf_dic['end_step'] = int(end_step)
    else:
      log_info.log_error('Input error: end_step should be integer, please check or reset analyze/rdf/end_step')
  else:
    rdf_dic['end_step'] = end_frame_id

  if ( 'r_increment' in rdf_dic.keys() ):
    r_increment = rdf_dic['r_increment']
    if ( data_op.eval_str(r_increment) == 2 ):
      rdf_dic['r_increment'] = float(r_increment)
    else:
      log_info.log_error('Input error: r_increment should be float, please check or reset analyze/rdf/r_increment')
  else:
    rdf_dic['r_increment'] = 0.1

  if ( 'box' in rdf_dic.keys() ):
    A_exist = 'A' in rdf_dic['box'].keys()
    B_exist = 'B' in rdf_dic['box'].keys()
    C_exist = 'C' in rdf_dic['box'].keys()
  else:
    log_info.log_error('Input error: no box, please set analyze/rdf/box')
    exit()

  if ( A_exist and B_exist and C_exist ):
    box_A = rdf_dic['box']['A']
    box_B = rdf_dic['box']['B']
    box_C = rdf_dic['box']['C']
  else:
    log_info.log_error('Input error: box setting error, please check analyze/rdf/box')
    exit()

  if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
    rdf_dic['box']['A'] = [float(x) for x in box_A]
  else:
    log_info.log_error('Input error: A vector of box wrong, please check analyze/rdf/box/A')
    exit()

  if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
    rdf_dic['box']['B'] = [float(x) for x in box_B]
  else:
    log_info.log_error('Input error: B vector of box wrong, please check analyze/rdf/box/B')
    exit()

  if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
    rdf_dic['box']['C'] = [float(x) for x in box_C]
  else:
    log_info.log_error('Input error: C vector of box wrong, please check analyze/rdf/box/C')
    exit()

  return rdf_dic

def check_spectrum_inp(spectrum_dic):

  '''
  check_spectrum_inp : check the input of spectrum.

  Args :
    spectrum_dic : dictionary
      spectrum_dic contains parameters for spectrum.
  Returns :
    spectrum_dic: dictionary
      spectrum_dic is the revised spectrum_dic.
  '''

  valid_type = ['general', 'water_mode', 'hydration_mode']
  if ( 'type' in spectrum_dic.keys() ):
    spec_type = spectrum_dic['type']
    if ( spec_type in valid_type ):
      pass
    else:
      log_info.log_error('Input error: %s is not supported, but check or reset analyze/power_spectrum/type' %(spec_type))
      exit()
  else:
    log_info('Input error: no type, please set analyze/power_spectrum/type')
    exit()

  if ( 'traj_vel_file' in spectrum_dic.keys() ):
    traj_vel_file = spectrum_dic['traj_vel_file']
    if ( os.path.exists(os.path.abspath(traj_vel_file)) ):
      spectrum_dic['traj_vel_file'] = os.path.abspath(traj_vel_file)
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(os.path.abspath(traj_vel_file), 'vel')
    else:
      log_info.log_error('Input error: %s file does not exist' %(traj_vel_file))
      exit()
  else:
    log_info.log_error('Input error: traj_vel_file, please set analyze/power_spectrum/traj_vel_file')
    exit()

  spec_type = spectrum_dic['type']
  if ( spec_type == 'water_mode' or spec_type == 'hydration_mode' ):
    if ( 'traj_coord_file' in spectrum_dic.keys() ):
      traj_coord_file = spectrum_dic['traj_coord_file']
      if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
        spectrum_dic['traj_coord_file'] = os.path.abspath(traj_coord_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
        exit()
    else:
      log_info.log_error('Input error: no traj_coord_file, please set analyze/power_spectrum/traj_coord_file')
      exit()

  if ( 'init_step' in spectrum_dic.keys() ):
    init_step = spectrum_dic['init_step']
    if ( data_op.eval_str(init_step) == 1 ):
      spectrum_dic['init_step'] = int(init_step)
    else:
      log_info.log_error('Input error: init_step should be integer, please check or reset analyze/power_spectrum/init_step')
      exit()
  else:
    spectrum_dic['init_step'] = start_frame_id

  if ( 'end_step' in spectrum_dic.keys() ):
    end_step = spectrum_dic['end_step']
    if ( data_op.eval_str(end_step) == 1 ):
      spectrum_dic['end_step'] = int(end_step)
    else:
      log_info.log_error('Input error: end_step should be integer, please check or reset analyze/power_spectrum/end_step')
      exit()
  else:
    spectrum_dic['end_step'] = end_frame_id

  if ( 'max_frame_corr' in spectrum_dic.keys() ):
    max_frame_corr = spectrum_dic['max_frame_corr']
    if ( data_op.eval_str(max_frame_corr) == 1 ):
      spectrum_dic['max_frame_corr'] = int(max_frame_corr)
    else:
      log_info.log_error('Input error: max_frame_corr should be integer, please check or reset analyze/power_spectrum/max_frame_corr')
      exit()
  else:
    spectrum_dic['max_frame_corr'] = int(frames_num/3)

  if ( 'start_wave' in spectrum_dic.keys() ):
    start_wave = spectrum_dic['start_wave']
    if ( data_op.eval_str(start_wave) == 1 or data_op.eval_str(start_wave) == 2 ):
      spectrum_dic['start_wave'] = float(start_wave)
    else:
      log_info.log_error('Input error: start_wave should be float, please check or reset analyze/power_spectrum/start_wave')
      exit()
  else:
    spectrum_dic['start_wave'] = 0

  if ( 'end_wave' in spectrum_dic.keys() ):
    end_wave = spectrum_dic['end_wave']
    if ( data_op.eval_str(end_wave) == 1 or data_op.eval_str(end_wave) == 2 ):
      spectrum_dic['end_wave'] = float(end_wave)
    else:
      log_info.log_error('Input error: end_wave should be float, please check or reset analyze/power_spectrum/end_wave')
      exit()
  else:
    spectrum_dic['end_wave'] = 0

  if ( 'normalize' in spectrum_dic.keys() ):
    normalize = spectrum_dic['normalize']
    if ( data_op.eval_str(normalize) == 1 ):
      if ( int(normalize) == 0 or int(normalize) == 1 ):
        pass
      else:
        log_info.log_error('Input error: normalize should be 0 or 1, please check or reset analyze/power_spectrum/normalize')
        exit()
    else:
      log_info.log_error('Input error: normalize should be 0 or 1, please check or reset analyze/power_spectrum/normalize')
      exit()
  else:
    spectrum_dic['normalize'] = 1

  if ( spec_type == 'general' or spec_type == 'water_mode' ):
    if ( 'atom_id' in spectrum_dic.keys() ):
      atom_id_list = data_op.get_id_list(spectrum_dic['atom_id'])
      spectrum_dic['atom_id'] = atom_id_list
    else:
      log_info.log_error('Input error: no atom_id, please set analyze/power_spectrum/atom_id')
      exit()
  else:
    if ( 'box' in spectrum_dic.keys() ):
      A_exist = 'A' in spectrum_dic['box'].keys()
      B_exist = 'B' in spectrum_dic['box'].keys()
      C_exist = 'C' in spectrum_dic['box'].keys()
    else:
      log_info.log_error('Input error: no box, please set analyze/power_spectrum/box')
      exit()

    if ( A_exist and B_exist and C_exist ):
      box_A = spectrum_dic['box']['A']
      box_B = spectrum_dic['box']['B']
      box_C = spectrum_dic['box']['C']
    else:
      log_info.log_error('Input error: box setting error, please check analyze/power_sepctrum/box')
      exit()

    if ( len(box_A) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_A) ):
      spectrum_dic['box']['A'] = [float(x) for x in box_A]
    else:
      log_info.log_error('Input error: A vector of box wrong, please check analyze/power_spectrum/box/A')
      exit()

    if ( len(box_B) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_B) ):
      spectrum_dic['box']['B'] = [float(x) for x in box_B]
    else:
      log_info.log_error('Input error: B vector of box wrong, please check analyze/power_spectrum/box/B')
      exit()

    if ( len(box_C) == 3 and all(data_op.eval_str(i) == 1 or data_op.eval_str(i) == 2 for i in box_C) ):
      spectrum_dic['box']['C'] = [float(x) for x in box_C]
    else:
      log_info.log_error('Input error: C vector of box wrong, please check analyze/power_spectrum/box/C')
      exit()

    if ( spec_type == 'hydration_mode' ):
      if ( 'hyd_shell_dist' in spectrum_dic.keys() ):
        hyd_shell_dist = spectrum_dic['hyd_shell_dist']
        if ( data_op.eval_str(hyd_shell_dist) == 1 or data_op.eval_str(hyd_shell_dist) == 2 ):
          spectrum_dic['hyd_shell_dist'] = float(hyd_shell_dist)
        else:
          log_info.log_error('Input error: hyd_shell_dist should be float, please check or reset analyze/power_spectrum/hyd_shell_dist')
          exit()
      else:
        log_info.log_error('Input error: no hyd_shell_dist, please set analyze/power_spectrum/hyd_shell_dist')
        exit()

      if ( 'dist_conv' in spectrum_dic.keys() ):
        dist_conv = spectrum_dic['dist_conv']
        if ( data_op.eval_str(dist_conv) == 1 or data_op.eval_str(dist_conv) == 2 ):
          spectrum_dic['dist_conv'] = float(dist_conv)
        else:
          log_info.log_error('Input error: dist_conv should be float, please check or reset analyze/power_spectrum/dist_conv')
          exit()
      else:
        spectrum_dic['dist_conv'] = 0.3

      if ( 'atom_type_pair' in spectrum_dic.keys() ):
        atom_type_pair = spectrum_dic['atom_type_pair']
        if ( len(atom_type_pair) == 2 and all(data_op.eval_str(x) == 0 for x in atom_type_pair) ):
          pass
        else:
          log_info.log_error('Input error: atom_type_pair should be 2 string, please check or reset analyze/power_spectrum/atom_type_pair')
          exit()
      else:
        log_info.log_error('Input error: no atom_type_pair, please set analyze/power_spectrum/atom_type_pair')
        exit()

  return spectrum_dic

def check_v_hartree_inp(v_hartree_dic):

  '''
  check_v_hartree_inp : check the input of v_hartree.

  Args :
    v_hartree_dic : dictionary
      v_hartree_dic contains parameters for v_hartree.
  Returns :
    v_hartree_dic: dictionary
      v_hartree_dic is the revised v_hartree_dic.
  '''

  v_hartree_dic = copy.deepcopy(v_hartree)

  if ( 'cube_file' in v_hartree_dic.keys() ):
    cube_file = v_hartree_dic['cube_file']
    if ( os.path.exists(os.path.abspath(cube_file)) ):
      v_hartree_dic['cube_file'] = os.path.abspath(cube_file)
    else:
      log_info.log_error('Input error: %s file does not exist' %(cube_file))
      exit()
  else:
    log_info.log_error('Input error: no charge cube file, please set analyze/v_hartree/cube_file')
    exit()

  if ( 'surface' in v_hartree_dic.keys() ):
    surface = v_hartree_dic['surface']
    if ( len(surface) == 3 and all(data_op.eval_str(x) == 1 for x in surface) ):
      v_hartree_dic['surface'] = [int(x) for x in surface]
    else:
      log_info.log_error('Input error: surface should be 3 integer, please check or set analyze/v_hartree/surface')
      exit()
  else:
    log_info.log_error('Input error: no surface, please set analyze/v_hartree/surface')
    exit()

  return v_hartree_dic

def check_arrange_data_inp(arrange_data_dic):

  '''
  check_arrange_data_inp : check the input of arrange_data.

  Args :
    arrange_data_dic : dictionary
      arrange_data_dic contains parameters for arrange_data.
  Returns :
    arrange_data_dic: dictionary
      arrange_data_dic is the revised arrange_data_dic.
  '''

  if ( 'temperature' in arrange_data_dic):
    temp_dic = arrange_data_dic['temperature']

    if ( 'traj_ener_file' in temp_dic.keys() ):
      traj_ener_file = temp_dic['traj_ener_file']
      if ( os.path.exists(os.path.abspath(traj_ener_file)) ):
        arrange_data_dic['temperature']['traj_ener_file'] = os.path.abspath(traj_ener_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_ener_file))
        exit()
    else:
      log_info.log_error('Input error: energy trajectory file, please set analyze/arranage_data/temperature/traj_ener_file')
      exit()

  #arrange potential energy
  elif ( 'potential' in arrange_data_dic ):
    pot_dic = arrange_data_dic['potential']

    if ( 'traj_ener_file' in pot_dic.keys() ):
      traj_ener_file = pot_dic['traj_ener_file']
      if ( os.path.exists(os.path.abspath(traj_ener_file)) ):
        arrange_data_dic['potential']['traj_ener_file'] = os.path.abspath(traj_ener_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_ener_file))
        exit()
    else:
      log_info.log_error('Input error: no energy trajectory file, please set analyze/arranage_data/potential/traj_ener_file')
      exit()

  #arrange mulliken charge
  elif ( 'mulliken' in arrange_data_dic ):
    mulliken_dic = arrange_data_dic['mulliken']

    if ( 'traj_mul_file' in mulliken_dic.keys() ):
      traj_mul_file = mulliken_dic['traj_mul_file']
      if ( os.path.exists(os.path.abspath(traj_mul_file)) ):
        arrange_data_dic['mulliken']['traj_mul_file'] = os.path.abspath(traj_mul_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_mul_file))
        exit()
    else:
      log_info.log_error('Input error: no mulliken trajectory file, please set analyze/arranage_data/mulliken/traj_mul_file')
      exit()

    if ( 'atom_id' in mulliken_dic.keys() ):
      arrange_data_dic['mulliken']['atom_id'] = data_op.get_id_list(mulliken_dic['atom_id'])
    else:
      log_info.log_error('Input error: no atom id, please set analyze/arranage_data/mulliken/atom_id')
      exit()

    if ( 'time_step' in mulliken_dic.keys() ):
      time_step = mulliken_dic['time_step']
      if ( data_op.eval_str(time_step) == 2 ):
        arrange_data_dic['mulliken']['time_step'] = float(time_step)
      else:
        log_info.log_error('Input error: time step should be float, please check or set analyze/arranage_data/mulliken/atom_id')
        exit()
    else:
      arrange_data_dic['mulliken']['time_step'] = 0.5

    if ( 'each' in mulliken_dic.keys() ):
      each = mulliken_dic['each']
      if ( data_op.eval_str(each) == 1 ):
        arrange_data_dic['mulliken']['each'] = int(each)
      else:
        log_info.log_error('Input error: each should be integer, please check or set analyze/arranage_data/mulliken/each')
        exit()
    else:
      arrange_data_dic['mulliken']['each'] = 1

  #arrange vertical energy
  elif ( 'vertical_energy' in arrange_data_dic ):
    vert_ene_dic = arrange_data_dic['vertical_energy']

    if ( 'traj_mix_ener_file' in vert_ene_dic.keys() ):
      traj_mix_ener_file = vert_ene_dic['traj_mix_ener_file']
      if ( os.path.exists(os.path.abspath(traj_mix_ener_file)) ):
        arrange_data_dic['vertical_energy']['traj_mix_ener_file'] = os.path.abspath(traj_mix_ener_file)
        blocks_num, base, pre_base, frames_num, each, start_id, end_id, time_step = \
        traj_info.get_traj_info(os.path.abspath(traj_mix_ener_file))
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_mix_ener_file))
        exit()
    else:
      log_info.log_error('Input error: no mix energy trajectory file, please set analyze/arranage_data/vertical_energy/traj_mix_ene_file')
      exit()

    if ( 'row_ox' in vert_ene_dic.keys() ):
      row_ox = vert_ene_dic['row_ox']
      if ( data_op.eval_str(row_ox) == 1 ):
        arrange_data_dic['vertical_energy']['row_ox'] = int(row_ox)
      else:
        log_info.log_error('Input error: row_ox should be integer, please check or reset analyze/arranage_data/vertical_energy/row_ox')
        exit()
    else:
      log_info.log_error('Input error: no row_ox, please set analyze/arranage_data/vertical_energy/row_ox')
      exit()

    if ( 'row_red' in vert_ene_dic.keys() ):
      row_red = vert_ene_dic['row_red']
      if ( data_op.eval_str(row_red) == 1 ):
        arrange_data_dic['vertical_energy']['row_red'] = int(row_red)
      else:
        log_info.log_error('Input error: row_red should be integer, please check or reset analyze/arranage_data/vertical_energy/row_red')
        exit()
    else:
      log_info.log_error('Input error: no row_red, please set analyze/arranage_data/vertical_energy/row_red')
      exit()

    if ( 'redox_type' in vert_ene_dic.keys() ):
      redox_type = vert_ene_dic['redox_type']
      if ( redox_type == 'oxidation' or redox_type == 'reduction' ):
        pass
      else:
        log_info.log_error('Input error: only oxidation and reduction are supported for redox_type')
        exit()
    else:
      log_info.log_error('Input error: no redox_type, please set analyze/arranage_data/vertical_energy/redox_type')
      exit()

    if ( 'slow_growth' in vert_ene_dic.keys() ):
      slow_growth = vert_ene_dic['slow_growth']
      if ( data_op.eval_str(slow_growth) == 1 ):
        if ( int(slow_growth) == 0 or int(slow_growth) == 1 ):
          arrange_data_dic['vertical_energy']['slow_growth'] = int(slow_growth)
        else:
          log_info.log_error('Input error: slow_growth should be 0 or 1, please check or reset analyze/arranage_data/vertical_energy/slow_growth')
          exit()
      else:
        log_info.log_error('Input error: slow_growth should be integer, please set analyze/arranage_data/vertical_energy/slow_growth')
        exit()
    else:
      arrange_data_dic['vertical_energy']['slow_growth'] = 0 #0 means no slow_growth

    slow_growth = arrange_data_dic['vertical_energy']['slow_growth']
    if ( slow_growth == 1 ):
      if ( 'increment' in vert_ene_dic.keys() ):
        increment = vert_ene_dic['increment']
        if ( data_op.eval_str(increment) == 2 ):
          arrange_data_dic['vertical_energy']['increment'] = float(increment)
        else:
          log_info.log_error('Input error: increment should be float, please check or reset analyze/arranage_data/vertical_energy/increment')
          exit()
      else:
        log_info.log_error('Input error: no increment, please set analyze/arranage_data/vertical_energy/increment')
        exit()

    if ( 'init_step' in vert_ene_dic.keys() ):
      init_step = vert_ene_dic['init_step']
      if ( data_op.eval_str(init_step) == 1 ):
        arrange_data_dic['vertical_energy']['init_step'] = int(init_step)
      else:
        log_info.log_error('Input error: init_step should be integer, please check or reset analyze/arranage_data/vertical_energy/init_step')
        exit()
    else:
      arrange_data_dic['vertical_energy']['init_step'] = start_id

    if ( 'end_step' in vert_ene_dic.keys() ):
      end_step = vert_ene_dic['end_step']
      if ( data_op.eval_str(end_step) == 1 ):
        arrange_data_dic['vertical_energy']['end_step'] = int(end_step)
      else:
        log_info.log_error('Input error: end_step should be integer, please check or reset analyze/arranage_data/vertical_energy/end_step')
        exit()
    else:
      arrange_data_dic['vertical_energy']['end_step'] = end_id

    if ( 'final_time_unit' in vert_ene_dic.keys() ):
      final_time_unit = vert_ene_dic['final_time_unit']
      if ( final_time_unit in ['fs', 'ps', 'ns'] ):
        pass
      else:
        log_info.log_error('Input error: final_time_unit wrong, please check or reset analyze/arranage_data/vertical_energy/final_time_unit')
        exit()
    else:
      arrange_data_dic['vertical_energy']['final_time_unit'] = 'fs'

  elif ( 'ti_force' in arrange_data_dic ):
    ti_force_dic = arrange_data_dic['ti_force']

    if ( 'traj_lag_file' in ti_force_dic.keys() ):
      traj_lag_file = ti_force_dic['traj_lag_file']
      if ( os.path.exists(os.path.abspath(traj_lag_file)) ):
        arrange_data_dic['ti_force']['traj_lag_file'] = os.path.abspath(traj_lag_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(traj_lag_file))
        exit()
    else:
      log_info.log_error('Input error: no traj_file, please set analyze/arrange_data/ti_force/traj_lag_file')
      exit()

    if ( 'stat_num' in ti_force_dic.keys() ):
      stat_num = ti_force_dic['stat_num']
      if ( data_op.eval_str(stat_num) == 1 ):
        arrange_data_dic['ti_force']['stat_num'] = int(stat_num)
      else:
        log_info.log_error('Input error: stat_num should be integer, please check or reset analyze/arrange_data/ti_force/stat_num')
        exit()
    else:
      arrange_data_dic['ti_force']['stat_num'] = 1

  return arrange_data_dic

def check_dp_test_inp(dp_test_dic):

  '''
  check_dp_test_inp : check the input of dp_test.

  Args :
    dp_test_dic : dictionary
      dp_test_dic contains parameters for dp_test.
  Returns :
    dp_test_dic: dictionary
      dp_test_dic is the revised dp_test_dic.
  '''

  if ( 'learn_type' in dp_test_dic.keys() ):
    learn_type = dp_test_dic['learn_type']
    if ( learn_type == 'supervised' or learn_type == 'active_learning' ):
      pass
    else:
      log_info.log_error('Input error: learn_type should be supervised or active_learning')
      exit()
  else:
    log_info.log_error('Input error: no learn_type, please set analyze/dp_test/learn_type')
    exit()

  learn_type = dp_test_dic['learn_type']

  if ( learn_type == 'supervised' ):
    if ( 'cp2k_frc_file' in dp_test_dic.keys() ):
      cp2k_frc_file = dp_test_dic['cp2k_frc_file']
      if ( os.path.exists(os.path.abspath(cp2k_frc_file)) ):
        dp_test_dic['cp2k_frc_file'] = os.path.abspath(cp2k_frc_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(cp2k_frc_file))
        exit()
    else:
      log_info.log_error('Input error: no cp2k_frc_file, please set analyze/dp_test/cp2k_frc_file')
      exit()

    if ( 'cp2k_pos_file' in dp_test_dic.keys() ):
      cp2k_pos_file = dp_test_dic['cp2k_pos_file']
      if ( os.path.exists(os.path.abspath(cp2k_pos_file)) ):
        dp_test_dic['cp2k_pos_file'] = os.path.abspath(cp2k_pos_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(cp2k_pos_file))
        exit()
    else:
      log_info.log_error('Input error: no cp2k_pos_file, please set analyze/dp_test/cp2k_pos_file')
      exit()

    if ( 'cp2k_cell_file' in dp_test_dic.keys() ):
      cp2k_cell_file = dp_test_dic['cp2k_cell_file']
      if ( os.path.exists(os.path.abspath(cp2k_cell_file)) ):
        dp_test_dic['cp2k_cell_file'] = os.path.abspath(cp2k_cell_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(cp2k_cell_file))
        exit()
    else:
      log_info.log_error('Input error: no cp2k_cell_file, please set analyze/dp_test/cp2k_cell_file')
      exit()

    if ( 'dpff_file' in dp_test_dic.keys() ):
      dpff_file = dp_test_dic['dpff_file']
      if ( os.path.exists(os.path.abspath(dpff_file)) ):
        dp_test_dic['dpff_file'] = os.path.abspath(dpff_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(dpff_file))
        exit()
    else:
      log_info.log_error('Input error: no dpff_file, please set analyze/dp_test/dpff_file')
      exit()

    if ( 'atom_label' in dp_test_dic.keys() ):
      atom_label = dp_test_dic['atom_label']
      atom_label_dic = OrderedDict()
      for i in range (len(atom_label)):
        lable_split = data_op.str_split(atom_label[i], ':')
        atom_label_dic[int(lable_split[0])] = lable_split[1]
      dp_test_dic['atom_lable'] = atom_lable_dic
    else:
      log_info.log_error('Input error: no atom_lable, please set analyze/dp_test/atom_lable')
      exit()

    supervised_test(cp2k_pos_file, cp2k_cell_file, cp2k_frc_file, dpff_file, atom_label_dic, work_dir)

  elif ( learn_type == 'active_learning' ):
    if ( 'lmp_traj_file' in dp_test_dic.keys() ):
      lmp_traj_file = dp_test_dic['lmp_traj_file']
      if ( os.path.exists(os.path.abspath(lmp_traj_file)) ):
        dp_test_dic['lmp_traj_file'] = os.path.abspath(lmp_traj_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(lmp_traj_file))
        exit()
    else:
      log_info.log_error('Input error: no lmp_traj_file, please set analyze/dp_test/lmp_traj_file')
      exit()

    if ( 'lmp_log_file' in dp_test_dic.keys() ):
      lmp_log_file = dp_test_dic['lmp_log_file']
      if ( os.path.exists(os.path.abspath(lmp_log_file)) ):
        dp_test_dic['lmp_log_file'] = os.path.abspath(lmp_log_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(lmp_log_file))
        exit()
    else:
      log_info.log_error('Input error: no lmp_log_file, please set analyze/dp_test/lmp_log_file')
      exit()

    if ( 'cp2k_inp_file' in dp_test_dic.keys() ):
      cp2k_inp_file = dp_test_dic['cp2k_inp_file']
      if ( os.path.exists(os.path.abspath(cp2k_inp_file)) ):
        dp_test_dic['cp2k_inp_file'] = os.path.abspath(cp2k_inp_file)
      else:
        log_info.log_error('Input error: %s file does not exist' %(cp2k_inp_file))
        exit()
    else:
      log_info.log_error('Input error: no cp2k_inp_file, please set analyze/dp_test/cp2k_inp_file')
      exit()

    if ( 'cp2k_exe' in dp_test_dic.keys() ):
      cp2k_exe = dp_test_dic['cp2k_exe']
      if ( os.path.exists(os.path.abspath(cp2k_exe)) ):
        dp_test_dic['cp2k_exe'] = os.path.abspath(cp2k_exe)
      else:
        log_info.log_error('Input error: cp2k executable file %s does not exist' %(cp2k_exe))
        exit()
    else:
      log_info.log_error('Input error: no cp2k executable file, please set analyze/dp_test/cp2k_exe')
      exit()

    if ( 'cp2k_mpi_num' in dp_test_dic.keys() ):
      cp2k_mpi_num = dp_test_dic['cp2k_mpi_num']
      if ( data_op.eval_str(cp2k_mpi_num) == 1 ):
        dp_test_dic['cp2k_mpi_num'] = int(cp2k_mpi_num)
      else:
        log_info.log_error('Input error: cp2k mpi numbers should be integer, please check or reset analyze/dp_test/cp2k_mpi_num')
        exit()
    else:
      dp_test_dic['cp2k_mpi_num'] = 1

  else:
    log_info.log_error('Input error: Could not recognize learn_type, please check or reset analyze/dp_test/learn_type')
    exit()

  if ( 'atom_label' in dp_test_dic.keys() ):
    atom_label = dp_test_dic['atom_label']
    atom_label_dic = OrderedDict()
    for i in range (len(atom_label)):
      label_split = data_op.str_split(atom_label[i], ':')
      atom_label_dic[int(label_split[0])] = label_split[1]
    dp_test_dic['atom_label'] = atom_label_dic
  else:
    log_info.log_error('Input error: no atom_label, please set analyze/dp_test/atom_label')
    exit()

  return dp_test_dic

def check_free_energy_inp(free_energy_dic):

  '''
  check_free_energy_inp : check the input of free_energy.

  Args :
    free_energy_dic : dictionary
      free_energy_dic contains parameters for free_energy.
  Returns :
    free_energy_dic: dictionary
      free_energy_dic is the revised free_energy_dic.
  '''

  free_energy_dic = copy.deepcopy(free_energy_dic)

  if ( 'method' in free_energy_dic.keys() ):
    method = free_energy_dic['method']
  else:
    log_info.log_error('Input error: no method, please set analyze/free_energy/method')
    exit()

  if ( method == 'ti' ):
    if ( 'ti_file' in free_energy_dic.keys() ):
      ti_file = free_energy_dic['ti_file']
      if ( os.path.exists(os.path.abspath(ti_file)) ):
        free_energy_dic['ti_file'] = os.path.abspath(ti_file)
      else:
        log_info.log_error('Input error: %s does not exist' %(ti_file))
        exit()
    else:
      log_info.log_error('Input error: no ti_file, please set analyze/free_energy/ti_file')
      exit()

  return free_energy_dic

def check_rmsd_inp(rmsd_dic):

  '''
  check_rmsd_inp : check the input of rmsd.

  Args :
    rmsd_dic : dictionary
      rmsd_dic contains parameters for rmsd.
  Returns :
    rmsd_dic: dictionary
      rmsd_dic is the revised rmsd_dic.
  '''

  if ( 'traj_coord_file' in rmsd_dic.keys() ):
    traj_coord_file = rmsd_dic['traj_coord_file']
    if ( os.path.exists(os.path.abspath(traj_coord_file)) ):
      rmsd_dic['traj_coord_file'] = os.path.abspath(traj_coord_file)
    else:
      log_info.log_error('Input error: %s file does not exist' %(traj_coord_file))
      exit()
  else:
    log_info.log_error('Input error: no coordination trajectory file, please set analyze/rmsd/traj_coord_file')
    exit()

  if ( 'atom_id' in rmsd_dic.keys() ):
    rmsd_dic['atom_id'] = data_op.get_id_list(rmsd_dic['atom_id'])
  else:
    log_info.log_error('Input error: no atom id, please set analyze/rmsd/atom_id')
    exit()

  if ( 'ref_frame' in rmsd_dic.keys() ):
    ref_frame = rmsd_dic['ref_frame']
    if ( data_op.eval_str(ref_frame) == 1 ):
      rmsd_dic['ref_frame'] = int(ref_frame)
    else:
      log_info.log_error('Input error: ref_frame should be integer, please check or reset analyze/rmsd/ref_frame')
      eixt()
  else:
    log_info.log_error('Input error: no reference frame, please set analyze/rmsd/ref_frame')
    exit()

  if ( 'compare_frame' in rmsd_dic.keys() ):
    compare_frame = rmsd_dic['compare_frame']
    rmsd_dic['compare_frame'] = data_op.get_id_list(compare_frame)
  else:
    log_info.log_error('Input error: no compare frame, please set analyze/rmsd/compare_frame')
    exit()

  return rmsd_dic

def check_time_correlation_inp(time_corr_dic):

  '''
  check_time_correlation_inp : check the input of time_correlation.

  Args :
    time_corr_dic : dictionary
      time_corr_dic contains parameters for time_corr.
  Returns :
    time_corr_dic: dictionary
      time_corr_dic is the revised time_corr_dic.
  '''

  if ( 'traj_file' in time_corr_dic.keys() ):
    traj_file = time_corr_dic['traj_file']
    if ( os.path.exists(os.path.abspath(traj_file)) ):
      time_corr_dic['traj_file'] = os.path.abspath(traj_file)
      atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
      traj_info.get_traj_info(os.path.abspath(traj_file), 'coord')
    else:
      log_info.log_error('Input error: %s file does not exists' %(traj_file))
      exit()
  else:
    log_info.log_error('Input error: no traj_file, please set analyze/time_correlation/traj_file')
    exit()

  if ( 'atom_id' in time_corr_dic.keys() ):
    atom_id = time_corr_dic['atom_id']
    time_corr_dic['atom_id'] = data_op.get_id_list(atom_id)
  else:
    log_info.log_error('Input error: no atom_id, please set analyze/time_correlation/atom_id')
    exit()

  if ( 'max_frame_corr' in time_corr_dic.keys() ):
    max_frame_corr = time_corr_dic['max_frame_corr']
    if ( data_op.eval_str(max_frame_corr) == 1 ):
      time_corr_dic['max_frame_corr'] = int(max_frame_corr)
    else:
      log_info.log_error('Input error: max_frame_corr should be integer, please check or reset analyze/time_correlation/max_frame_corr')
      exit()
  else:
    time_corr_dic['max_frame_corr'] = int(frames_num/3)

  if ( 'init_step' in time_corr_dic.keys() ):
    init_step = time_corr_dic['init_step']
    if ( data_op.eval_str(init_step) == 1 ):
      time_corr_dic['init_step'] = int(init_step)
    else:
      log_info.log_error('Input error: init_step should be integer, please check or reset analyze/time_correlation/init_step')
      exit()
  else:
    time_corr_dic['init_step'] = start_frame_id

  if ( 'end_step' in time_corr_dic.keys() ):
    end_step = time_corr_dic['end_step']
    if ( data_op.eval_str(end_step) == 1 ):
      time_corr_dic['end_step'] = int(end_step)
    else:
      log_info.log_error('Input error: end_step should be integer, please check or reset analyze/time_correlation/end_step')
      exit()
  else:
    time_corr_dic['end_step'] = end_frame_id

  if ( 'normalize' in time_corr_dic.keys() ):
    normalize = time_corr_dic['normalize']
    if ( data_op.eval_str(normalize) == 1 ):
      if ( int(normalize) == 0 or int(normalize) == 1 ):
        time_corr_dic['normalize'] = int(init_step)
      else:
        log_info.log_error('Input error: only 0 and 1 are supported for normalize, please check or reset analyze/time_correlation/normalize')
        exit()
    else:
      log_info.log_error('Input error: init_step should be integer, please check or reset analyze/time_correlation/normalize')
      exit()
  else:
    time_corr_dic['normalize'] = 1

  return time_corr_dic
