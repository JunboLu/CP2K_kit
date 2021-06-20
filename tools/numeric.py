#!/usr/bin/env python

import numpy as np
from math import factorial

def get_float_power(f_temp):

  '''
  get_float_power : get the power value of a float in a scientific notation.

  Args :
    f_temp : float
      Example : 0.001
  Returns :
    int, the power value of a float in a scientific notation.
    Example : -3
  '''

  f_temp_str = '{:e}'.format(f_temp)
  f_temp_str_split = f_temp_str.split('e')

  return int(f_temp_str_split[1])

def get_as_num_string(f_temp):

  '''
  get_as_num_string : transform a float with 8 decimal into string.

  Args :
    f_temp : float
      Example : 0.001
  Returns :
    y : string
      Example : '0.00100000'
  '''

  y = '{:.8f}'.format(f_temp)

  return y

def get_cosine(a_arr, b_arr):

  '''
  get_cosine : get cosine similarity between two data sets

  Args :
    a_arr : 1-d float array
    b_arr : 1-d float array
  Return :
    cos_value : float
      cos_value is cosine value between two vector.
  '''

  a_dot = np.dot(a_arr,a_arr)
  b_dot = np.dot(b_arr,b_arr)
  a_dot_sqrt = np.sqrt(a_var)
  b_dot_sqrt = np.sqrt(b_var)
  ab_dot = np.dot(a_arr,b_arr)

  cos_value = ab_dot/(a_dot_sqrt*b_dot_sqrt)

  return cos_value

def get_euclid_dist(a_arr, b_arr):

  '''
  get_euclid_dist : get euclid distance between two data sets

  Args :
    a_arr : 1-d float array
    b_arr : 1-d float array
  Return :
    euclid_dist : float
      euclid_dist is euclid distance between two vector.
  '''

  sum_value = 0.0
  for i in range(len(a_arr)):
    sum_value = sum_value + (a_arr[i]-b_arr[i])**2

  euclid_dist = np.sqrt(sum_value/(len(a_arr)))

  return euclid_dist

def get_abs_list(list_a, list_b):

  '''
  get_abs_list : get abs of difference between list_a and list_b

  Args :
    list_a : 1-d float list
    list_b : 1-d float list
  '''

  list_c = []
  for i in range(len(list_a)):
    list_c.append(abs(list_a[i]-list_b[i]))

  return list_c

def get_corr_coeff(a_arr, b_arr):

  '''
  get_corr_coeff : get correlation coefficient between two data sets

  Args :
    a_arr : 1-d float array
    b_arr : 1-d float array
  Return :
    corr_value : float
      corr_value is correlation coefficient between two data sets.
  '''

  #Get mean value
  a_mean = a_arr.mean()
  b_mean = b_arr.mean()

  a_de_mean = [ a_i - a_mean for a_i in a_arr ]
  b_de_mean = [ b_i - b_mean for b_i in b_arr ]
  a_de_mean_arr = np.array(a_de_mean)
  b_de_mean_arr = np.array(b_de_mean)

  a_var = np.dot(a_de_mean_arr, a_de_mean_arr)/(len(a_de_mean_arr)-1)
  b_var = np.dot(b_de_mean_arr, b_de_mean_arr)/(len(b_de_mean_arr)-1)

  a_devi = np.sqrt(a_var)
  b_devi = np.sqrt(b_var)

  ab_cov = np.dot(a_de_mean_arr, b_de_mean_arr)/(len(a_de_mean_arr)-1)

  if ( a_devi > 0.0 and b_devi > 0.0 ):
    corr_value = ab_cov/(a_devi*b_devi)

  return corr_value

def savitzky_golay(y, window_size, order, deriv=0, rate=1):

  '''
  savitzky_golay : fit the curve by using savitzky-golay method.
  '''

  window_size = np.abs(np.int(window_size))
  order = np.abs(np.int(order))
  order_range = range(order+1)
  half_window = (window_size -1) // 2

  b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
  m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)

  firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
  lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
  y = np.concatenate((firstvals, y, lastvals))

  return np.convolve( m[::-1], y, mode='valid')

if __name__ == '__main__':
  from CP2K_kit.tools import numeric
  a = np.array([100.0,5.0,3.0,0.8,1.0,0.9,0.2,0.4,0.6,0.7,0.3,0.5,0.1,0.2,0.3,0.4])
  b = np.array([1000.0,10.0,8.0,0.9,0.9,1.0])
  c = np.array([99.0,4.0,2.9,0.01,0.02,0.03,-1,-2,-3,-4,-5,-6,-10,-15,-11,-12])
  d = np.array([1,0.9,0.9,1,1,1,1,1,1,1])
  e = np.array([10,9.9,9.9,10,10,10,10,10,10,10])
  f = np.array([0.1,0.2])
  corr_coeff = numeric.get_corr_coeff(e, d)
#  corr_coeff = numeric.get_euclid_dist(f, f[::-1])
  print (corr_coeff)

