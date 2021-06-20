#! /usr/env/bin python

import sys
from CP2K_kit.deepff import active

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])

active.kernel(work_dir, inp_file)
