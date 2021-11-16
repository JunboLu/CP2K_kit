#! /bin/env/python

import linecache
import os

gau_num=3
cp2k_version='cp2k-6.1'
pre = ['0.1','0.2','0.3','0.4','0.5','0.6']
conv_1 =  30.0
conv_2 =  1.0

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
      if (gau_num == 3):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[2])<20.0 and (n[2]-n[1])>0.2 and (n[1]-n[0])>0.1):
          print (i,(j+1),(k+1),n[0],n[1],n[2])
      elif (gau_num == 4):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[3])<20.0 and (n[3]-n[2])>0.2 and (n[2]-n[1])>0.1 and (n[1]-n[0])>0.05):
          print (i,(j+1),(k+1),n[0],n[1],n[2],n[3])
      elif (gau_num == 5):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[4])<20.0 and (n[4]-n[3])>0.2 and (n[3]-n[2])>0.1 and (n[2]-n[1])>0.05):
          print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4])
      elif (gau_num == 6):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[5])<20.0 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.01):
          print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5])
      elif (gau_num == 7):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[6])<conv_1 and float(n[6])>conv_2 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.02 and (n[1]-n[0])>0.015):
          print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6])
      elif (gau_num ==8):
        if (len(n) == gau_num and float(n[0])<0.1 and float(n[7])<20.0 and (n[7]-n[6])>0.3 and (n[6]-n[5])>0.3 and (n[5]-n[4])>0.2 and (n[4]-n[3])>0.1 and (n[3]-n[2])>0.05 and (n[2]-n[1])>0.01):
          print (i,(j+1),(k+1),n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7])

