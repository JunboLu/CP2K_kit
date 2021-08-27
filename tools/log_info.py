#! /usr/env/bin python

from CP2K_kit.tools import data_op

def log_logo():

  '''
  log_log: print cp2k_kit logo
  '''

  logo = '''
********************************************************************************
********************************************************************************
**                                                                            **
**                                                                            **
**       ######    #####    ######  ##   ##       ##   ##  ##  ########       **
**     ##        ##     ##      ##  ##  ##        ##  ##   ##     ##          **
**     ##        ##     ##      ##  ## ##         ## ##    ##     ##          **
**     ##        #######    ######  ####   ###### ####     ##     ##          **
**     ##        ##         ##      ## ##         ## ##    ##     ##          **
**       ######  ##         ##      ##  ##        ##  ##   ##     ##          **
**               ##         ######  ##   ##       ##   ##  ##     ##          **
**                                                                            **
**                        CP2K-KIT : A helper of CP2K                         **
**                        Email : lujunbo15@gmail.com                         **
**                                                                            **
********************************************************************************
'''

  print (logo, flush=True)

def log_traj_info(atoms_num, frames_num, each, start_frame_id, end_frame_id, time_step):

  '''
  log_traj_info: print trajectory information

  Args:
    atoms_num: int
      atoms_num is the number of atoms in the system.
    frames_num: int
      frames_num is the number of frames in trajectory file.
    each: int
      each is printing frequency of the trajectory.
    start_frame_id: int
      start_frame_id is the starting frame id in trajectory file.
    end_frame_id: int
      end_frame_id is the endding frame id in trajectory file.
    time_step: float
      time_step is time step of md. Its unit is fs in CP2K_kit.
  Returns:
    none
  '''

  print ('Trajectory info'.center(80,'*'), flush=True)
  print ('The number of atoms for the system is %d' %(atoms_num), flush=True)
  print ('The number of frames in the trajectory is %d' %(frames_num), flush=True)
  print ('The printing frequency of the trajectory is %d' %(each), flush=True)
  print ('The trajectory is started from %d step' %(start_frame_id), flush=True)
  print ('The trajectory is ended in %d step' %(end_frame_id), flush=True)
  print ('The time step is %f fs\n' %(time_step), flush=True)

def log_error(error, error_lable='Error'):

  '''
  log_error: print error information

  Args:
    error: string
      error is the error information
    error_lable: string
      error_lable is the lable of error
  Returns:
    none
  '''

  print (error_lable.center(80,'*'), flush=True)
  print (data_op.str_wrap(error,80), flush=True)
  print ('\n', flush=True)
