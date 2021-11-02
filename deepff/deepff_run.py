#! /usr/env/bin python

import sys
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.tools import read_input
from CP2K_kit.deepff import active
from CP2K_kit.deepff import dp_test

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
deepff_type = str(sys.argv[3])

log_info.log_logo()

print (data_op.str_wrap('DEEPFF| PROGRAM STARTED IN %s' %(work_dir), 80), flush=True)
print ('DEEPFF| INPUT FILE NAME %s\n' %(inp_file), flush=True)

if ( deepff_type == 'active_model_devi' or deepff_type == 'active_dp_test' ):
  active.kernel(work_dir, inp_file, deepff_type)
elif ( deepff_type == 'dp_test' ):
  job_type_param = read_input.dump_info(work_dir, inp_file, [deepff_type])
  dp_test.dp_test_run(job_type_param[0], work_dir)
