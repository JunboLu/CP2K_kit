#! /usr/env/bin python

import sys
from CP2K_kit.tools import read_input
from CP2K_kit.model import droplet

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
model_type = str(sys.argv[3])

model_type_param = read_input.dump_info(work_dir, inp_file, [model_type])

if ( model_type == 'droplet' ):
  droplet.droplet_run(model_type_param[0], work_dir)

