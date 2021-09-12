#! /usr/env/bin python

import sys
from CP2K_kit.tools import log_info
from CP2K_kit.tools import data_op
from CP2K_kit.deepff import active

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])

log_info.log_logo()

print (data_op.str_wrap('DEEPFF| PROGRAM STARTED IN %s' %(work_dir), 80), flush=True)
print ('DEEPFF| INPUT FILE NAME %s\n' %(inp_file), flush=True)

active.kernel(work_dir, inp_file)
