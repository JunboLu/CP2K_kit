#! /usr/env/bin python

import copy
from CP2K_kit.tools import data_op

def gen_atom_inp(work_dir, gth_pp_opt_param, weight_1):

  '''
  gen_atom_inp : generate atom input file used to optimize gth pp.

  Args :
    work_dir : string
      work_dir is working directory of CP2K_kit.
    gth_pp_opt_param : dictionary
      gth_pp_opt_param contains keywords used to generate atom input file.
    weight_1: 1-d float list
      weight_1 is the initial weight.
  Returns :
    none
  '''

  element = gth_pp_opt_param['element']
  elec_config = gth_pp_opt_param['elec_config']
  elec_core_config = gth_pp_opt_param['elec_core_config']
  all_elec_method = gth_pp_opt_param['all_elec_method']
  xc_func = gth_pp_opt_param['xc_func']
  relat_method = gth_pp_opt_param['relat_method']

  inp_file = ''.join((work_dir, '/atom.inp'))
  atom_file = open(inp_file, 'w')

  atom_file.write('&GLOBAL\n')
  atom_file.write('  PROGRAM_NAME ATOM\n')
  atom_file.write('&END GLOBAL\n')

  atom_file.write('&ATOM\n')

  element_line = ''.join(('  ELEMENT ', element, '\n'))
  atom_file.write(element_line)
  atom_file.write('  RUN_TYPE PSEUDOPOTENTIAL_OPTIMIZATION\n')

  config = ''
  for i in range(len(elec_config)):
    config = ' '.join((config, elec_config[i]))

  elec_config_line = ''.join(('  ELECTRON_CONFIGURATION', config, '\n'))
  atom_file.write(elec_config_line)

  if ( isinstance(elec_core_config, list) ):
    core_config = ''
    for i in range(len(elec_core_config)):
      core_config = ' '.join((core_config, elec_core_config[i]))
    core_config_line = ''.join(('  CORE', core_config, '\n'))
  else:
    core_config_line = ''.join(('  CORE  ', elec_core_config, '\n'))
  atom_file.write(core_config_line)

  val_config = copy.deepcopy(elec_config)
  for i in elec_config:
   if ( i in elec_core_config ):
     val_config.remove(i)

  val_orbit = []
  val_orbit_occ = []
  for i in range(len(val_config)):
    val_orbit.append(val_config[i][1])
    val_orbit_occ.append(int(val_config[i][2:len(val_config[i])]))

  val_orbit_type = data_op.list_replicate(val_orbit)

  if ( 'f' in val_orbit_type ):
    max_ang = '3'
  else:
    if ( 'd' in val_orbit_type ):
      max_ang = '2'
    else:
      if ( 'p' in val_orbit_type ):
        max_ang = '1'
      else:
        if ( 's' in val_orbit_type ):
          max_ang = '0'

  max_ang_line = ''.join(('  MAX_ANGULAR_MOMENTUM ', max_ang, '\n'))
  atom_file.write(max_ang_line)

  atom_file.write('  &METHOD\n')

  method_type_line = ''.join(('    METHOD_TYPE  ', all_elec_method.upper(), '\n'))
  atom_file.write(method_type_line)

  relat_method_line = ''.join(('    RELATIVISTIC ', relat_method.upper(), '\n'))
  atom_file.write(relat_method_line)

  if ( all_elec_method == 'kohn-sham' ):
    atom_file.write('    &XC\n')

    xc_func_line = ''.join(('      &XC_FUNCTIONAL ', xc_func.upper(), '\n'))
    atom_file.write(xc_func_line)

    atom_file.write('      &END XC_FUNCTIONAL\n')
    atom_file.write('    &END XC\n')

  atom_file.write('  &END METHOD\n')
  atom_file.write('  &OPTIMIZATION\n')
  atom_file.write('    EPS_SCF 1.e-7\n')
  atom_file.write('  &END OPTIMIZATION\n')
  atom_file.write('  &AE_BASIS\n')
  atom_file.write('    BASIS_TYPE GEOMETRICAL_GTO\n')
  atom_file.write('  &END AE_BASIS\n')
  atom_file.write('  &PP_BASIS\n')
  atom_file.write('     BASIS_TYPE GEOMETRICAL_GTO\n')
  atom_file.write('  &END PP_BASIS\n')
  atom_file.write('  &POTENTIAL\n')
  atom_file.write('    PSEUDO_TYPE GTH\n')
  atom_file.write('    POTENTIAL_FILE_NAME  ../GTH-PARAMETER\n')

  if ( all_elec_method == 'kohn-sham' ):
    pot_name = ''.join(('GTH-', xc_func.upper(), '-q', str(sum(val_orbit_occ))))
  else:
    pot_name = ''.join(('GTH-', all_elec_method, '-q', str(sum(val_orbit_occ))))

  pot_name_line = ''.join(('    POTENTIAL_NAME ', pot_name, '\n'))
  atom_file.write(pot_name_line)

  atom_file.write('    CONFINEMENT_TYPE  BARRIER\n')
  atom_file.write('    CONFINEMENT 200. 4.0 12.0\n')
  atom_file.write('  &END POTENTIAL\n')

  atom_file.write('  &POWELL\n')
  atom_file.write('    ACCURACY   1.e-14\n')
  atom_file.write('    STEP_SIZE  0.01\n')
  atom_file.write('    MAX_INIT   200\n')
  atom_file.write('    MAX_FUN    50\n')
  atom_file.write('    STEP_SIZE_SCALING  0.90\n')
  atom_file.write('    WEIGHT_PSIR0 10.0\n')
  atom_file.write('    TARGET_POT_SEMICORE      [eV]      0.003000\n')
  atom_file.write('    TARGET_POT_VALENCE       [eV]      0.000300\n')
  atom_file.write('    TARGET_POT_VIRTUAL       [eV]      0.003000\n')
  atom_file.write('    WEIGHT_POT_NODE                    10.0\n')
  atom_file.write('    WEIGHT_POT_SEMICORE                %f\n' %(weight_1[0]))
  atom_file.write('    WEIGHT_POT_VALENCE                 %f\n' %(weight_1[1]))
  atom_file.write('    WEIGHT_POT_VIRTUAL                 %f\n' %(weight_1[2]))
  atom_file.write('    SEMICORE_LEVEL       [eV]         15.0\n')
  atom_file.write('  &END\n')
  atom_file.write('&END ATOM\n')

  atom_file.close()

  if ( all_elec_method == 'kohn-sham' ):
    return element, sum(val_orbit_occ), xc_func.upper()
  else:
    return element, sum(val_orbit_occ), all_elec_method.upper()
