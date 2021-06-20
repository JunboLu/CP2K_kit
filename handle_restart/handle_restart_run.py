#! /usr/env/bin python

import sys
from CP2K_kit.tools import read_input
from CP2K_kit.handle_restart import handle_restart_traj

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
run_type = str(sys.argv[3])

job_type_param = read_input.dump_info(work_dir, inp_file, [run_type])
print (job_type_param)

handle_restart_traj.handle_restart_run(job_type_param[0])
