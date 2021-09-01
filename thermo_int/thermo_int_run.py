#! /usr/env/bin python

import sys
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info
from CP2K_kit.tools import read_input
from CP2K_kit.thermo_int import mix_force_eval
from CP2K_kit.thermo_int import constraint_md

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
thermo_type = str(sys.argv[3])

log_info.log_logo()

print (data_op.str_wrap('THERMO_INT| PROGRAM STARTED IN %s' %(work_dir), 80), flush=True)
print ('THERMO_INT| Input file name %s\n' %(inp_file), flush=True)

thermo_type_param = read_input.dump_info(work_dir, inp_file, [thermo_type])

if ( thermo_type == 'mix_force_eval' ):
  mix_force_eval.mix_force_eval_run(thermo_type_param[0], work_dir)
elif ( thermo_type == 'constraint_md' ):
  constraint_md.constraint_md_run(thermo_type_param[0], work_dir)
