#!/usr/bin/env python

def get_atom_mass(atom_name):

  '''
  get_atom_mass : get atom number and mass.

  Args :
    atom_name : string or 1-d string list
      atom_name is atom name. Example : 'O'/['O','H']
  Returns :
    atom_num : int or 1-d int list
      atom_num is atom number. Example : 6/[6,1]
    atom_mass : float or 1-d float list
      atom_mass is atom mass. Example : 16.00/[16.00,1.01]
  '''

  atom_id_mass = {'H':(1,1.01), 'He':(2,4.00), 'Li':(3,6.94), 'Be':(4,9.01), 'B':(5,10.81),
               'C':(6,12.01), 'N':(7,14.01), 'O':(8,16.00), 'F':(9,19.00), 'Ne':(10,20.18),
               'Na':(11,22.99), 'Mg':(12,24.31), 'Al':(13,26.98), 'Si':(14,28.09), 'P':(15,30.97),
               'S':(16,32.07), 'Cl':(17,35.45), 'Ar':(18,39.95), 'K':(19,39.10), 'Ca':(20,40.08),
               'Sc':(21,44.96), 'Ti':(22,47.87), 'V':(23,50.94), 'Cr':(24,52.00), 'Mn':(25,54.94),
               'Fe':(26,55.85), 'Co':(27,58.93), 'Ni':(28,58.69), 'Cu':(29,63.55), 'Zn':(30,65.38),
               'Ga':(31,69.72), 'Ge':(32,72.64), 'As':(33,74.92), 'Se':(34,78.96), 'Br':(35,79.90),
               'Kr':(36,83.80), 'Rb':(37,85.47), 'Sr':(38,87.62), 'Y':(39,88.91), 'Zr':(40,91.22),
               'Nb':(41,92.91), 'Mo':(42,95.96), 'Tc':(43,98), 'Ru':(44,101.07), 'Rh':(45,102.91),
               'Pd':(46,106.42), 'Ag':(47,107.87), 'Cd':(48,112.41), 'In':(49,114.82), 'Sn':(50,118.71),
               'Sb':(51,121.76), 'Te':(52,127.60), 'I':(53,126.90), 'Xe':(54,131.29), 'Cs':(55,132.91),
               'Ba':(56,137.33), 'La':(57,138.91), 'Ce':(58,140.12), 'Pr':(59,140.91), 'Nd':(60,144.24),
               'Pm':(61,145), 'Sm':(62,150.36), 'Eu':(63,151.96), 'Gd':(64,157.25), 'Tb':(65,158.93),
               'Dy':(66,162.50), 'Ho':(67,164.93), 'Er':(68,167.26), 'Tm':(69,168.93), 'Yb':(70,173.05),
               'Lu':(71,174.97), 'Hf':(72,178.49), 'Ta':(73,180.95), 'W':(74,183.84), 'Re':(75,186.21),
               'Os':(76,190.23), 'Ir':(77,192.22), 'Pt':(78,1985.08), 'Au':(79,196.97), 'Hg':(80,200.59),
               'Tl':(81,204.38), 'Pb':(82,207.2), 'Bi':(83,208.98), 'Po':(84,209), 'At':(85,210),
               'Rn':(86,222), 'Fr':(87,223), 'Ra':(88,226), 'Ac':(89,227), 'Th':(90,232.04),
               'Pa':(91,231.04), 'U':(92,238.03), 'Np':(93,237), 'Pu':(94,244), 'Am':(95,243),
               'Cm':(96,247), 'Bk':(97,247), 'Cf':(98,251), 'Es':(99,252), 'Fm':(100,257),
               'Md':(101,258), 'No':(102,259), 'Lr':(103,262)}

  #atom could be a string list, and we will return float list.
  if ( isinstance(atom_name, list) ):
    atom_num = []
    atom_mass = []

    for i in range(len(atom_name)):
      if ( atom_name[i] in atom_id_mass ):
        atom_num.append(atom_id_mass[atom_name[i]][0])
        atom_mass.append(atom_id_mass[atom_name[i]][1])
      else:
        print ("Doesn't find %s atom in atom_id_mass" % (atom_name[i]))
        exit()

  #atom could be a string, and we will return float.
  elif ( isinstance(atom_name, str) ):
    if ( atom_name in atom_id_mass ):
      atom_num = atom_id_mass[atom_name][0]
      atom_mass = atom_id_mass[atom_name][1]
    else:
      print ("Doesn't find %s atom in atom_id_mass" %(atom_name))
      exit()

  return atom_num, atom_mass

def get_atom_cov_radius(atom):

  '''
  get_atom_cov_radius : get atom covalence radius.

  Args :
    atom : string or 1-d string list
      atom is atom name. Example : 'O'/['O','H']
  Returns :
    atom_cov_rad : float or 1-d float list
      atom_cov_rad is atom covalence radius. Example : 0.66/[0.66, 0.3]
  '''

  #Unit of covalence radius is Angstrom.
  atom_id_cov_radius = {'H':0.3, 'Li':1.23, 'Be':0.89, 'B':0.81,
               'C':0.77, 'N':0.7, 'O':0.66, 'F':0.64,
               'Na':1.57, 'Mg':1.36, 'Al':1.25, 'Si':1.17, 'P':1.10,
               'S':1.04, 'Cl':0.99, 'K':2.03, 'Ca':1.74,
               'Sc':1.44, 'Ti':1.32, 'V':1.22, 'Cr':1.18, 'Mn':1.18,
               'Fe':1.16, 'Co':1.16, 'Ni':1.15, 'Cu':1.17, 'Zn':1.25,
               'Ga':1.25, 'Ge':1.22, 'As':1.21, 'Se':1.17, 'Br':1.14,
               'Rb':2.16, 'Sr':1.91, 'Y':1.62, 'Zr':1.45,
               'Nb':1.34, 'Mo':1.30, 'Tc':1.27, 'Ru':1.25, 'Rh':1.25,
               'Pd':1.28, 'Ag':1.34, 'Cd':1.41, 'In':1.50, 'Sn':1.40,
               'Sb':1.41, 'Te':1.37, 'I':1.33, 'Cs':2.53,
               'Ba':1.98, 'La':1.69, 'Ce':1.65, 'Pr':1.65, 'Nd':1.64,
               'Pm':1.63, 'Sm':1.64, 'Eu':1.85, 'Gd':1.61, 'Tb':1.59,
               'Dy':1.59, 'Ho':1.58, 'Er':1.57, 'Tm':1.56, 'Yb':1.74,
               'Lu':1.56, 'Hf':1.44, 'Ta':1.34, 'W':1.30, 'Re':1.28,
               'Os':1.26, 'Ir':1.27, 'Pt':1.30, 'Au':1.34, 'Hg':1.44,
               'Tl':1.55, 'Pb':1.54, 'Bi':1.46, 'Po':1.46, 'At':1.46,
               'Ac':1.65, 'Th':1.65,'Pa':1.65, 'U':1.42, 'Np':1.42, 'Pu':1.42, 'Am':1.42,
               'Cm':1.42, 'Bk':1.42, 'Cf':1.42, 'Es':1.42, 'Fm':1.42,
               'Md':1.42, 'No':1.42, 'Lr':1.42}

  #atom could be string list.
  if ( isinstance(atom, list) ):
    atom_cov_rad = []

    for i in range(len(atom)):
      if ( atom[i] in atom_id_cov_radius ):
        atom_cov_rad.append(atom_id_cov_radius[atom[i]])
      else:
        print ("Doesn't find %s atom in atom_id_cov_radius" % (atom[i]))
        exit()

  #atom could be string.
  elif ( isinstance(atom, str) ):
    if ( atom in atom_id_cov_radius ):
      atom_cov_rad = atom_id_cov_radius[atom]
    else:
      print ("Doesn't find %s atom in atom_id_cov_radius" %(atom))
      exit()

  return atom_cov_rad

def get_q_info(atom, q_type='large'):

  '''
  get_q_info : get the number of valence electrons for a atom.
    This function is used to be in generating CP2K input file.

  Args :
    atom : string
      atom is atom name. Example : 'O'
    q_type : string
      q_type is the core type for GTH pseudopotential.
      There are three choices: 'large', 'medium', 'small'
  Returns :
    q_num : string
      q_num is the number of valence electrons. Example : 'q6'
  '''

  atom_q_large = {'H':'q1', 'He':'q2', 'Li':'q3', 'Be':'q4', 'B':'q3', 'C':'q4',
                  'N':'q5', 'O':'q6', 'F':'q7', 'Ne':'q8', 'Na':'q1', 'Mg':'q2',
                  'Al':'q3', 'Si':'q4', 'P':'q5', 'S':'q6', 'Cl':'q7', 'Ar':'q8',
                  'K':'q1', 'Ca':'q2', 'Sc':'q3', 'Ti':'q4', 'V':'q5', 'Cr':'q6',
                  'Mn':'q7', 'Cu':'q11', 'Zn':'q12', 'Ga':'q3', 'Ge':'q4', 'As':'q5',
                  'Br':'q7', 'Kr':'q8', 'Ru':'q8', 'Rh':'q9', 'Pd':'q10', 'Ag':'q11',
                  'Cd':'q12', 'In':'q3', 'Sn':'q4', 'Sb':'q5', 'Te':'q6', 'I':'q7',
                  'Xe':'q8', 'Cs':'q1', 'Ba':'q2', 'La':'q11', 'Ce':'q12', 'Pr':'q13',
                  'Nd':'q14', 'Pm':'q15', 'Sm':'q16', 'Eu':'q17', 'Gd':'q18', 'Ta':'q5',
                  'W':'q6', 'Re':'q7', 'Os':'q8', 'Ir':'q9', 'Pt':'q10', 'Au':'q11', 'Hg':'q12',
                  'Tl':'q3', 'Pb':'q4', 'Bi':'q5', 'Po':'q6', 'At':'q7', 'Rn':'q8', 'Ac':'q11',
                  'Th':'q12', 'Pa':'q13', 'U':'q14', 'Np':'q15', 'Pu':'q16', 'Am':'q17', 'Cm':'q18',
                  'Bk':'q19', 'Cf':'q20', 'Es':'q21', 'Fm':'q22', 'Md':'q23', 'No':'q24', 'Lr':'q25' }

  atom_q_medium = {'Na':'q9', 'Mg':'q10', 'K':'q9', 'Ca':'q10', 'Sc':'q11', 'Ti':'q12',
                   'V':'q13', 'Cr':'q14', 'Mn':'q15', 'Fe':'q16', 'Co':'q17', 'Ni':'q18',
                   'Cu':'q19', 'Zn':'q20', 'Ca':'q13', 'Rb':'q9', 'Sr':'q10', 'Y':'q11',
                   'Zr':'q12', 'Nb':'q13', 'Mo':'q14', 'Tc':'q15', 'Ru':'q16', 'Rh':'q17',
                   'Pd':'q18', 'Ag':'q19', 'In':'q13', 'Cs':'q9', 'Ba':'q10', 'Tb':'q29',
                   'Dy':'q30', 'Ho':'q31', 'Er':'q32', 'Tm':'q33', 'Yb':'q34', 'Lu':'q35',
                   'Hf':'q12', 'Ta':'q13', 'W':'q14', 'Re':'q15', 'Os':'q16', 'Ir':'q17',
                   'Pt':'q18', 'Au':'q19', 'Tl':'q13', 'Pb':'q14', 'Bi':'q15', 'Ac':'q21',
                   'Th':'q22', 'Pa':'q23', 'U':'q24', 'Np':'q25', 'Pu':'q26', 'Am':'q27',
                   'Cm':'q28', 'Bk':'q29', 'Cf':'q30', 'Es':'q31', 'Fm':'q32', 'Md':'q33',
                   'No':'q34', 'Lr':'q35'}

  if ( q_type == 'large'):
    if ( atom in atom_q_large ):
      q_num = atom_q_large[atom]
    else:
      print ("Doesn't find %s atom in atom_q_large" % (atom))
      exit()

  elif ( q_type == 'medium'):
    if ( atom in atom_q_medium ):
      q_num = atom_q_medium[atom]
    else:
      #If we cannot find atom in atom_q_medium, we will search it in atom_q_medium
      if ( atom in atom_q_large ):
        q_num = atom_q_large[atom]
      else:
        print ("Doesn't find %s atom in atom_q_large or atom_q_medium" % (atom))
        exit()

  return atom_q_large[atom]

if __name__ == '__main__':
  from CP2K_kit.tools import atom
  atom_num, atom_mass = atom.get_atom_mass(['O', 'H'])
  print (atom_num, atom_mass)
  atom_cov_rad = atom_rev.get_atom_cov_radius(['O','H'])
  print (atom_cov_rad)
