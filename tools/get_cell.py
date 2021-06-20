#! /usr/env/bin python

import numpy as np

def get_cell_const(a_vec, b_vec, c_vec):

  '''
  get_cell_const : get cell constant from cell vectors.

  Args :
    a_vec : 1-d float array, dim = 3
      a_vec is the cell vector of a.
      Example : array([12.42, 0.0, 0.0])
    b_vec : 1-d float array, dim = 3
      b_vec is the cell vector of b.
      Example : array([0.0, 12.42, 0.0])
    c_vec : 1-d float array, dim = 3
      c_vec is the cell vector of c
      Example : array([0.0, 0.0, 12.42])
  Returns :
    a : float
      a is the length of a_vec.
    b : float
      b is the length of b_vec.
    c : float
      c is the length of c_vec.
    alpha : float
      alpha is the angle between b_vec and c_vec.
    beta : float
      beta is the angle between a_vec and c_vec.
    gamma : float
      gamma is the angle between a_vec and b_vec.
  '''

  a = np.sqrt(np.dot(a_vec,a_vec))
  b = np.sqrt(np.dot(b_vec,b_vec))
  c = np.sqrt(np.dot(c_vec,c_vec))

  ab = np.dot(a_vec,b_vec)
  ac = np.dot(a_vec,c_vec)
  bc = np.dot(b_vec,c_vec)

  #numpy returns rad
  alpha = np.arccos(bc/(b*c))
  beta = np.arccos(ac/(a*c))
  gamma = np.arccos(ab/(a*b))

  return a, b, c, alpha, beta, gamma

def get_triclinic_cell(a_vec, b_vec, c_vec):

  '''
  get_triclinic_cell : get triclinic cell from conventional cell vectors.

  Args :
    a_vec : 1-d float array, dim = 3
      a_vec is the cell vector of a.
      Example : array([12.42, 0.0, 0.0])
    b_vec : 1-d float array, dim = 3
      b_vec is the cell vector of b.
      Example : array([0.0, 12.42, 0.0])
    c_vec : 1-d float array, dim = 3
      c_vec is the cell vector of c.
      Example : array([0.0, 0.0, 12.42])
  Returns :
    tri_a_vec : 1-d float array, dim = 3
      tri_a_vec is the cell vector of a in triclinic cell.
      Example : array([Lx, 0.0, 0.0])
    tri_b_vec : 1-d float array, dim = 3
      tri_b_vec is the cell vector of b in triclinic cell.
      Example : array([xy, Ly, 0.0])
    tri_c_vec : 1-d float array, dim = 3
      Example : array([xz, yz, Lz])
  '''

  a, b, c, alpha, beta, gamma = get_cell_const(a_vec, b_vec, c_vec)

  tri_a_vec = np.zeros(len(a_vec))
  tri_b_vec = np.zeros(len(b_vec))
  tri_c_vec = np.zeros(len(c_vec))

  tri_a_vec[0] = a
  tri_b_vec[0] = b*np.cos(gamma)
  tri_b_vec[1] = b*np.sin(gamma)
  tri_c_vec[0] = c*np.cos(beta)
  tri_c_vec[1] = c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
  temp = 1.0+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.cos(alpha)**2-np.cos(beta)**2-np.cos(gamma)**2
  tri_c_vec[2] = c*np.sqrt(temp)/np.sin(gamma)

  return tri_a_vec, tri_b_vec, tri_c_vec

def get_triclinic_cell_six(box_param):

  '''
  get_triclinic_cell_six : get triclinic cell from six parameters
  Args :
    box_param : 1-d float list, dim = 6
      Example : [Lx, Ly, Lz, xy, xz, yz]
  Returns :
    tri_a_vec : 1-d float array, dim = 3
      Example : array([Lx, 0.0, 0.0])
    tri_b_vec : 1-d float array, dim = 3
      Example : array([xy, Ly, 0.0])
    tri_c_vec : 1-d float array, dim = 3
      Example : array([xz, yz, Lz])
  '''

  tri_a_vec = np.zeros(3)
  tri_b_vec = np.zeros(3)
  tri_c_vec = np.zeros(3)

  tri_a_vec[0] = box_param[0]
  tri_b_vec[0] = box_param[3]
  tri_b_vec[1] = box_param[1]
  tri_c_vec[0] = box_param[4]
  tri_c_vec[1] = box_param[5]
  tri_c_vec[2] = box_param[2]

  return tri_a_vec, tri_b_vec, tri_c_vec

if __name__ == '__main__':
  from CP2K_kit.tools import get_triclinic
  a_vec = np.array([1.0,0.0,0.0])
  b_vec = np.array([1.0,1.0,0.0])
  c_vec = np.array([0.0,0.0,1.0])
#  a, b, c, alpha, beta, gamma = get_triclinic.get_cell_const(a_vec, b_vec, c_vec)
#  print (a, b, c, alpha, beta, gamma)
  new_a_vec, new_b_vec, new_c_vec = get_cell.get_triclinic_cell(a_vec, b_vec, c_vec)
  print (new_a_vec, new_b_vec, new_c_vec)

