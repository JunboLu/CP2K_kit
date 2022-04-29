#! /usr/env/bin python

import linecache

force_ti_file = open('force_ti', 'w')
for i in range(20000):
  line = linecache.getline('Lu.LagrangeMultLog', i*2+1)
  line_split = line.split(' ')
  force_ti = line_split[len(line_split)-1].strip('\n')
  distance = 1.03 + 0.000025*i
  force_ti_file.write('%f %s\n' %(distance, force_ti))

force_ti_file.close()
