#! /usr/env/bin python

from CP2K_kit.thermo_int import redox_pka

work_dir = str(sys.argv[1])
inp_file = str(sys.argv[2])
thermo_type = str(sys.argv[3])

thermo_type_param = read_input.dump_info(work_dir, inp_file, [thermo_type])

if ( thermo_type == 'redox_pka' ):
  redox_pka.run_redox_pka(thermo_type_param, work_dir)


