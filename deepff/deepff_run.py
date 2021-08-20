#! /usr/env/bin python

import sys
from CP2K_kit.tools import log_info
from CP2K_kit.tools import list_dic_op
from CP2K_kit.deepff import active

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])

log_info.log_logo()

print (list_dic_op.str_wrap('DEEPFF| PROGRAM STARTED IN %s' %(work_dir), 80), flush=True)
print ('DEEPFF| Input file name %s\n' %(inp_file), flush=True)

active.kernel(work_dir, inp_file)
