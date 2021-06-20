#! /usr/env/bin python

import sys
from CP2K_kit.tools import traj_info
from CP2K_kit.tools import read_input
from CP2K_kit.analyze import rdf
from CP2K_kit.analyze import rmsd
from CP2K_kit.analyze import center
from CP2K_kit.analyze import geometry
from CP2K_kit.analyze import diffusion
from CP2K_kit.analyze import spectrum
from CP2K_kit.analyze import free_energy
from CP2K_kit.analyze import arrange_data
from CP2K_kit.analyze import time_correlation

#We add a new keyword: analyze_job. We will use this keyword to assign job. 
work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
analyze_job = str(sys.argv[3])

job_type_param = read_input.dump_info(work_dir, inp_file, [analyze_job])

#People could do traj_info at first as it will help them to know the 
#basis information of the trajectory.
if ( analyze_job == 'traj_info' ):
  traj_info_param = job_type_param[0]
  atoms_num, base, pre_base, frames_num, each, start_frame_id, end_frame_id, time_step = \
  traj_info.get_traj_info(traj_info_param['traj_file'])
  print ('The number of atoms is %d' % (atoms_num))
  print ('The number of frames is %d' % (frames_num))
  print ('The printing frequency is %d' % (each))
  print ('The staring frame is %d' % (start_frame_id))
  print ('The endding frame is %d' % (end_frame_id))
  print ('The time step is %f fs' % (time_step))

elif ( analyze_job == 'center' ):
  center.center_run(job_type_param[0], work_dir)

elif ( analyze_job == 'geometry' ) :
  geometry.geometry_run(job_type_param[0], work_dir)

elif ( analyze_job == 'arrange_data' ):
  arrange_data.arrange_data_run(job_type_param[0], work_dir)

elif ( analyze_job == 'diffusion' ):
  diffusion.diffusion_run(job_type_param[0], work_dir)

elif ( analyze_job == 'rdf' ):
  rdf.rdf_run(job_type_param[0], work_dir)

elif ( analyze_job == 'rmsd' ):
  rmsd.rmsd_run(job_type_param[0], work_dir)

elif ( analyze_job == 'time_correlation' ):
  time_correlation.time_corr_run(job_type_param[0], work_dir)

elif ( analyze_job == 'power_spectrum' ):
  spectrum.power_spectrum_run(job_type_param[0], work_dir)
