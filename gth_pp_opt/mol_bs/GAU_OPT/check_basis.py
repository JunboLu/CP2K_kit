#! /env/bin/python

import os
import sys
import copy
import linecache

###########################################################
#!!!Change parameters here!!!
gau_num=6
cp2k_version='cp2k-6.1'
orb_type = 'GTO'
conv_1 = 20.0
conv_2 = 4.0
###########################################################

pre = ['0.1','0.2','0.3','0.4','0.5','0.6']
small_gau_valid = [0.1, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.055, \
                   0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01]
small_gau_stat_0 = []
small_gau_stat_1 = []
small_gau_stat_2 = []
small_gau_stat_3 = []
small_gau_stat_4 = []
small_gau_stat_5 = []
small_gau_stat_6 = []
small_gau_stat_7 = []
small_gau_stat_8 = []
small_gau_stat_9 = []
small_gau_stat_10 = []
small_gau_stat_11 = []
small_gau_stat_12 = []
small_gau_stat_13 = []
small_gau_stat_14 = []
small_gau_stat_15 = []
small_gau_stat_16 = []
small_gau_stat_17 = []
small_gau_stat_18 = []
large_gau_stat_0 = []
large_gau_stat_1 = []
large_gau_stat_2 = []
large_gau_stat_3 = []
large_gau_stat_4 = []
large_gau_stat_5 = []
large_gau_stat_6 = []
large_gau_stat_7 = []
large_gau_stat_8 = []
large_gau_stat_9 = []
large_gau_stat_10 = []
large_gau_stat_11 = []
large_gau_stat_12 = []
large_gau_stat_13 = []
large_gau_stat_14 = []
large_gau_stat_15 = []
large_gau_stat_16 = []
large_gau_stat_17 = []
large_gau_stat_18 = []

def judge_gau(small_gau_stat, large_gau_stat, n):

  state = False
  num = len(n)
  small_gau_stat_copy = copy.deepcopy(small_gau_stat)
  large_gau_stat_copy = copy.deepcopy(large_gau_stat)
  if ( len(small_gau_stat) == 0 ) :
    if ( num == 3 and (n[2]-n[1])>0.2 and (n[1]-n[0])>0.05 ):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
    elif ( num == 4 and (n[3]-n[2])>0.2 and (n[2]-n[1])>0.1 and (n[1]-n[0])>0.05 ):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
    elif ( num == 5 and (n[4]-n[3])>0.2 and (n[3]-n[2])>0.1 and (n[2]-n[1])>0.05 and (n[1]-n[0])>0.02 ):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
    elif ( num == 6 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 \
           and (n[1]-n[0])>0.01 ):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
    elif ( num == 7 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and \
         (n[2]-n[1])>0.02 and (n[1]-n[0])>0.01 ):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
    elif ( num == 8 and (n[7]-n[6])>0.3 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and \
         (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 and (n[1]-n[0])>0.01):
      state = True
      small_gau_stat_copy.append(n[0])
      large_gau_stat_copy.append(n[num-1])
  else:
    small_gau_minus_n_0 = [abs(x-n[0]) for x in small_gau_stat]
    large_gau_minus_n_num = [abs(x-n[num-1]) for x in large_gau_stat]
    if ( min(large_gau_minus_n_num) >= 0.1 ):
      if ( num == 3 and (n[2]-n[1])>0.2 and (n[1]-n[0])>0.05 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])
      elif ( num == 4 and (n[3]-n[2])>0.2 and (n[2]-n[1])>0.1 and (n[1]-n[0])>0.05 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])
      elif ( num == 5 and (n[4]-n[3])>0.2 and (n[3]-n[2])>0.1 and (n[2]-n[1])>0.05 and (n[1]-n[0])>0.02 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])
      elif ( num == 6 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 \
             and (n[1]-n[0])>0.01 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])
      elif ( num == 7 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and \
           (n[2]-n[1])>0.02 and (n[1]-n[0])>0.01 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])
      elif ( num == 8 and (n[7]-n[6])>0.3 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and \
           (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 and (n[1]-n[0])>0.01 ):
        state = True
        small_gau_stat_copy.append(n[0])
        large_gau_stat_copy.append(n[num-1])

  return state, small_gau_stat_copy, large_gau_stat_copy

def log_out(state, gau_num, n, i, j, k):

  if ( state and gau_num == 5 ):
    print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4], flush=True)
  elif ( state and gau_num == 6 ):
    print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5], flush=True)
  elif ( state and gau_num == 7 ):
    print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6], flush=True)
  elif ( state and gau_num == 8 ):
    print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7], flush=True)

for i in pre:
  for j in range(40):
    for k in range(25):
      opt_basis_file = 'pre_'+i+'/'+'range_'+str(j+1)+'/'+'exp_'+str(k+1)+'/'+'OPT_BASIS'
      if (os.path.exists(opt_basis_file)):
        if ( cp2k_version == 'cp2k-4.1' ):
          b = linecache.getline(opt_basis_file, 4)
          str_cont = '*' in b
          if (not str_cont == True):
            c = []
            n = []
            d = b.split(' ')
            for m in range(len(d)):
              if (d[m] != ''):
                c.append(d[m])
            c[len(c)-1] = c[len(c)-1].strip('\n')
            if (len(c) == gau_num+1):
              for z in range(len(c)-1):
                n.append(float(c[z+1]))
            n.sort()
        else:
          n = []
          for l in range(gau_num):
            if ( orb_type == 'GTO' ):
              b = linecache.getline(opt_basis_file, l+2)
            elif ( orb_type == 'GEO_GTO' ):
              b = linecache.getline(opt_basis_file, l+4)
            d = b.split(' ')
            c = []
            for m in range(len(d)):
              if (d[m] != ''):
                c.append(d[m])
            if ( l == 0 ):
              n.append(float(c[1]))
            else:
              n.append(float(c[0]))
          n.sort()
        if ( len(n) == gau_num and float(n[gau_num-1]) < conv_1 and float(n[gau_num-1]) >= conv_2 ):
          stat_index = 100
          for index, value in enumerate(small_gau_valid):
            if ( abs(float(n[0])-value) <= 0.001 ):
              stat_index = index
              break
          if ( stat_index == 0 ):
            state, small_gau_stat_0, large_gau_stat_0 = judge_gau(small_gau_stat_0, large_gau_stat_0, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 1 ):
            state, small_gau_stat_1, large_gau_stat_1 = judge_gau(small_gau_stat_1, large_gau_stat_1, n)  
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 2 ):
            state, small_gau_stat_2, large_gau_stat_2 = judge_gau(small_gau_stat_2, large_gau_stat_2, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 3 ):
            state, small_gau_stat_3, large_gau_stat_3 = judge_gau(small_gau_stat_3, large_gau_stat_3, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 4 ):
            state, small_gau_stat_4, large_gau_stat_4 = judge_gau(small_gau_stat_4, large_gau_stat_4, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 5 ):
            state, small_gau_stat_5, large_gau_stat_5 = judge_gau(small_gau_stat_5, large_gau_stat_5, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 6 ):
            state, small_gau_stat_6, large_gau_stat_6 = judge_gau(small_gau_stat_6, large_gau_stat_6, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 7 ):
            state, small_gau_stat_7, large_gau_stat_7 = judge_gau(small_gau_stat_7, large_gau_stat_7, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 8 ):
            state, small_gau_stat_8, large_gau_stat_8 = judge_gau(small_gau_stat_8, large_gau_stat_8, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 9 ):
            state, small_gau_stat_9, large_gau_stat_9 = judge_gau(small_gau_stat_9, large_gau_stat_9, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 10 ):
            state, small_gau_stat_10, large_gau_stat_10 = judge_gau(small_gau_stat_10, large_gau_stat_10, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 11 ):
            state, small_gau_stat_11, large_gau_stat_11 = judge_gau(small_gau_stat_11, large_gau_stat_11, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 12 ):
            state, small_gau_stat_12, large_gau_stat_12 = judge_gau(small_gau_stat_12, large_gau_stat_12, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 13 ):
            state, small_gau_stat_13, large_gau_stat_13 = judge_gau(small_gau_stat_13, large_gau_stat_13, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 14 ):
            state, small_gau_stat_14, large_gau_stat_14 = judge_gau(small_gau_stat_14, large_gau_stat_14, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 15 ):
            state, small_gau_stat_15, large_gau_stat_15 = judge_gau(small_gau_stat_15, large_gau_stat_15, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 16 ):
            state, small_gau_stat_16, large_gau_stat_16 = judge_gau(small_gau_stat_16, large_gau_stat_16, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 17 ):
            state, small_gau_stat_17, large_gau_stat_17 = judge_gau(small_gau_stat_17, large_gau_stat_17, n)
            log_out(state, gau_num, n, i, j, k)
          elif ( stat_index == 18 ):
            state, small_gau_stat_18, large_gau_stat_18 = judge_gau(small_gau_stat_18, large_gau_stat_18, n)
            log_out(state, gau_num, n, i, j, k)
