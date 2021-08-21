#! /usr/bin/env python

import re
import copy
import numpy as np

#This module defines several operation functions to list, integer, dictionary and string.
#These functions are heavilly used.

#Functions for list
def gen_list(start, end, incre):

  '''
  gen_list : generate a list from start, end and increment.

  Args :
    start : int
      start is starting value of the list.
    end : int
      end is the endding value of the list.
    incre : int
      incre is the increment of the list
  Returns :
    list_num : 1-d int list
      list_num is the final generated list.
  '''

  num = int((end-start)/incre)+1
  list_num = []
  for i in range(num):
    list_num.append(start+i*incre)

  return list_num

def comb_list_2_str(list_temp, space_char):

  '''
  comb_list_2_str : combine list element as string
  Transfer [1,2,3,4,5] to '1 2 3 4 5'

  Args :
    list_temp : 1-d list
      list_temp is the list needed to be combined.
    space_char : string
      space_char is the space element.
  Returns :
    comb_str : string
      comb_str is the combined string.
  '''

  comb_str = ''
  for i in range(len(list_temp)):
    comb_str = comb_str + str(list_temp[i]) + space_char

  return comb_str.strip(space_char)

def list_reshape(list_temp):

  #Here the dummy list_temp is 2-d list, for the higher dim
  #list, the systematical method is needed.

  '''
  list_reshape : reshape the given list

  Args :
    list_temp : 2-d list
    list_true : 1-d list
  '''

  list_true = []
  for i in range(len(list_temp)):
    for j in range(len(list_temp[i])):
      list_true.append(list_temp[i][j])

  return list_true

def order_list(list_temp, order_index):

  #The dims of list_temp and order_index should be same.

  '''
  order_list : order the list with an ordered index

  Args :
    list_temp : 1-d list
    order_index : 1-d int list
  Returns :
    list_true : 1-d list
  '''

  list_true = []
  for i in range(len(order_index)):
    list_true.append(list_temp[order_index[i]])

  return list_true

def list_num_stat(list_temp, value_cut, compare):

  '''
  list_num_stat : get the number of elements in a list which are greater or less than value_cut

  Args :
    list_temp : 1-d float list
    value_cut : float
      value_cut is the cutoff value.
    compare : string
      compare is the comparing type. Two choices : 'greater' and 'less'
  Returns :
    num : int
      num is the number of elements in a list which are greater or less than value_cut.
  '''

  num = 0
  if ( compare == 'greater' ):
    for i in range(len(list_temp)):
      if ( list_temp[i] >= value_cut ):
        num = num + 1
  if ( compare == 'less' ):
    for i in range(len(list_temp)):
      if ( list_temp[i] <= value_cut ):
        num = num + 1

  return num

def list_split(list_temp, n):

  '''
  list_split : divide the list into several parts, and each parts have n elements at most.

  Args :
    list_temp : 1-d list.
      Example : [1, 2, 3, 4, 5, 6, 7, 8, 9]
    n : int
      n is the number of elements in the sub-list. Example : 5
  Return :
    sub-lists. Example : [1, 2, 3, 4, 5], [6, 7, 8, 9]
  '''

  for i in range(0, len(list_temp), n):
    yield list_temp[i:i+n]

def list_replicate(list_temp, return_num=False):

  '''
  list_replicate : remove the same element in a list.

  Args :
    list_temp : 1-d list
      Example : ['a', 'b', 'c', 'a', 'b']
    return_num : bool
      return_num is whether we need return the number of elements in a list.
  Returns :
    list_type : 1-d list
      list_type is the final list where the same elements are removed.
      Example : ['a', 'b', 'c']
    num : 1-d int list
      num is the number of elenments in a list.
      Example : [2, 2, 1]
  '''

  list_type = []
  for x in list_temp:
    if x not in list_type:
      list_type.append(x)

  num = []
  for i in range(len(list_type)):
    num.append(0)

  for i in range(len(list_type)):
    for j in range(len(list_temp)):
      if ( list_temp[j] == list_type[i] ):
        num[i] = num[i] + 1

  if return_num == True:
    return list_type, num
  elif return_num == False:
    return list_type

def list_order(list_temp, order, return_index=False):

  '''
  list_order : order a list.

  Args :
    list_temp : 1-d int or float list
    order : string
      order is the order type. Two choices : 'ascend' and 'descend'.
    return_index : bool
      return_index is whether we need to return index.
  Returns:
    list_order : 1-d int or float list
      list_order is the final ordered list.
    list_order_index : 1-d int list
      list_order_index is the index of element in ordered list.
  '''

  if ( order == 'ascend' ):
    list_order = np.sort(np.array(list_temp))
    list_order_index = np.argsort(np.array(list_temp))
  if ( order == 'descend' ):
    list_order = np.sort(np.array(list_temp))[::-1]
    list_order_index = np.argsort(-np.array(list_temp))

  #list_order and list_order_index are array, we return list.
  if ( return_index == False ):
    return list(list_order)
  else:
    return list(list_order), list(list_order_index)

#Functions for string
def eval_str(str_tmp):

  '''
  eval_str: evaluate string

  Args:
    str_tmp: string
      str_tmp is the string needed to be evaluated.
  Returns:
    if the string is evaluated as string, return 0.
    if the string is evaluated as int, return 1.
    if the string is evaluated as float, return 2.
  '''

  try:
    int(str_tmp)
    return 1
  except ValueError:
    try:
      float(str_tmp)
      return 2
    except ValueError:
      return 0

def str_wrap(str_tmp, max_width, indent_char=''):

  '''
  str_wrap : wrap the string

  Args :
    str_temp : string
      str_temp is the string needed to be wrapped.
    max_width : int
      max_width is the maximum width.
  Returns :
    str_wrap : string
  '''

  if ( len(str_tmp) <= max_width ):
    return indent_char+str_tmp
  else:
    max_width = max_width - len(indent_char)
    list_tmp = [str_tmp[i:i+max_width] for i in range(0, len(str_tmp), max_width)]
    for i in range(len(list_tmp)):
      list_tmp[i] = indent_char + list_tmp[i]
    str_wrap = '\n'.join(list_tmp)
    return str_wrap

def str_split(str_temp, space_char):

  '''
  str_split : split the string and then return a list.

  Args :
    str_temp : string
      str_temp is the string needed to be splitted.
      Example : 'a b c'
    space_char : string
      space_char is the splitting element.
      Example : ' '
  Returns :
    list_temp : 1-d string list
      list_temp is the final splited list.
      Example : ['a', 'b', 'c']
  '''

  str_temp_split = str_temp.split(space_char)
  list_temp = []
  for i in range(len(str_temp_split)):
    if ( str_temp_split[i] != ''):
      list_temp.append(str_temp_split[i])

  return list_temp

def get_str_num(str_temp):

  '''
  get_str_num : get the number string in a string.

  Args :
    str_temp : string
      Example : '123/q'
  Return :
    string
      Example : '123'
  '''

  return re.sub('\D', '', str_temp)

def get_id_list(id_str):

  '''
  get_id_list : get a list from a string

  Args :
    id_str : string or list
      Example : "1-10", "1 8 10 30" "1-3 8 13 15-18"
      Example : ['1','2'], ['1','1-10']
  Returns :
    id_list : 1-d int list
      Example : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
                [1, 8, 10, 30];
                [1, 2, 3, 8, 13, 15, 16, 17, 18]
  '''

  if ( isinstance(id_str, list) ):
    a = copy.deepcopy(id_str)
  elif ( isinstance(id_str, str) ):
    a = str_split(id_str, ' ')

  id_list = []
  for i in range(len(a)):
    if ( '-' in a[i] and ':' not in a[i] ):
      b = a[i].index('-')
      start_id = int(a[i][:b])
      end_id = int(a[i][(b+1):])
      for j in range(end_id-start_id+1):
        id_list.append(start_id+j)
    elif ( '-' in a[i] and ':' in a[i] ):
      b = a[i].index('-')
      c = a[i].index(':')
      start_id = int(a[i][:b])
      end_id = int(a[i][b+1:c])
      increment = int(a[i][c+1:])
      num = int((end_id-start_id)/increment)+1
      for j in range(num):
        id_list.append(start_id+j*increment)
    else:
      id_list.append(int(a[i]))

  return id_list

def str_to_bool(str_tmp):

  #This function is easy to understand, so we do not prepare examples.

  '''
  str_to_bool : transfrom string variable to bool variable.

  Args :
    str_temp : string
  Returns :
    bool_str : bool
  '''

  if ( str_tmp.upper() == 'TRUE' ):
    bool_str = True
    return bool_str
  elif ( str_tmp.upper() == 'FALSE' ):
    bool_str = False
    return bool_str
  else:
    return str_tmp

#Functions for dictionary
def get_dic_keys(dic, value):

  #In general, one could get value for the given key in a dictionary. 
  #But we hope to get key from value in this function.
  #Please note that different keys will have the same value, so be careful
  #when using this function.

  '''
  get_dic_keys : get the keys from the value

  Args :
    dic : dictionary
      Example : {'a':1, 'b':2, 'c':3}
    value : int, float, string, list or dictionary
      Example : 1
  Returns:
    Example : 'a'
  '''

  return [k for k,v in dic.items() if v == value]

#Functions for integer:
def int_split(int_tmp, n):

  '''
  int_split : split one integer as n parts
  '''

  assert n > 0
  quotient = int(int_tmp/n)

  remainder = int_tmp%n
  if ( remainder > 0 ):
    return [quotient]*(n-remainder) + [quotient+1]*remainder
  if ( remainder < 0 ):
    return [quotient-1]*(-remainder) + [quotient]*(n+remainder)

  return [quotient]*n

if __name__ == '__main__':
  from CP2K_kit.tools import list_dic_op
  #Test comb_list_2_str:
  list_temp = [1,2,3,4,5]
  comb_str = list_dic_op.comb_list_2_str(list_temp, ' ')
  print (comb_str)
