#! /bin/env/python

import linecache
import os

pre = ['0.1','0.2','0.3','0.4','0.5','0.6']
gau_num = 7
conv_1 =  30.0
conv_2 =  1.0

for i in pre:
  for j in range(40):
    for k in range(25):
      input = 'pre_'+i+'/'+'range_'+str(j+1)+'/'+'exp_'+str(k+1)+'/'+'OPT_BASIS'
      if (os.path.exists(input)):
        b = linecache.getline(input, 4)
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
      if (gau_num == 6): 
        if (len(c) == gau_num+1 and float(n[0])<0.1 and float(n[5])<20.0 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.01):
          print i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5]
      elif (gau_num == 7):
        if (len(c) == gau_num+1 and float(n[0])<0.1 and float(n[6])<conv_1 and float(n[6])>conv_2 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 and (n[1]-n[0])>0.015):
          print i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6]
      elif (gau_num ==8):
        if (len(c) == gau_num+1 and float(n[0])<0.1 and float(n[7])<20.0 and (n[7]-n[6])>0.3 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.01):
          print i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]

