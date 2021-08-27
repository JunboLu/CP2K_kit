#! /usr/env/bin python

import subprocess
import linecache
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op

def read_gpuinfo(work_dir, gpuinfo_file):

  '''
  read_gpuinfo: read gpu information from gpu info file.

  Args:
    work_dir: string
      work_dir is the working directory.
    gpuinfo_file: string
      gpuinfo_file is a file containing gpu information.
  Returns:
    device: 1-d string list
      device is the gpu device name.
    usage: 1-d float list
      usage is the memory use of gpu device.
  '''

  device = []
  usage = []

  cmd_a = "grep -n %s %s" % ("'|===============================+'", gpuinfo_file)
  a = call.call_returns_shell(work_dir, cmd_a)
  if ( a != [] ):
    a_int = int(a[0].split(':')[0])
    cmd_b = "grep %s %s" % ("'+-------------------------------+'", gpuinfo_file)
    b = call.call_returns_shell(work_dir, cmd_b)
    device_num = len(b)
    for i in range(device_num):
      line_1 = linecache.getline(gpuinfo_file, a_int+i*4+1)
      line_1_split = data_op.str_split(line_1, ' ')
      device.append(line_1_split[1])
      line_2 = linecache.getline(gpuinfo_file, a_int+i*4+2)
      line_2_split = data_op.str_split(line_2, ' ')
      mem_used = float(line_2_split[8][0:line_2_split[8].index('M')])
      mem_tot = float(line_2_split[10][0:line_2_split[10].index('M')])
      usage.append(mem_used/mem_tot)

  return device, usage

def analyze_gpu(host, ssh, work_dir):

  '''
  analyze_gpu: analyze the gpu node information

  Args:
    host: 1-d string list
      host is the computational nodes.
    ssh: bool
      ssh is whether we need to ssh.
    work_dir: string
      work_dir is the working directory.
  Returns:
    device: 2-d string list
      device is the gpu device name.
    usage: 2-d float list
      usage is the memory use of gpu device.
  '''

  device = []
  usage = []

  for i in range(len(host)):
    gpuinfo_file = ''.join((work_dir, '/gpuinfo_', host[i]))
    if ssh:
      check_gpu = '''
#! /bin/bash

ssh %s
nvidia-smi > %s
exit
''' %(host[i], gpuinfo_file)

    else:
      check_gpu = '''
#! /bin/bash

nvidia-smi > %s
''' %(gpuinfo_file)

    check_gpu_file = ''.join((work_dir, '/check_gpu.sh'))
    with open(check_gpu_file, 'w') as f:
      f.write(check_gpu)

    subprocess.run('chmod +x check_gpu.sh', cwd=work_dir, shell=True)
    subprocess.run("bash -c './check_gpu.sh'", cwd=work_dir, shell=True)

    device_i, usage_i = read_gpuinfo(work_dir, gpuinfo_file)
    device.append(device_i)
    usage.append(usage_i)

  cmd_1 = 'rm check_gpu.sh'
  call.call_simple_shell(work_dir, cmd_1)
  cmd_2 = 'rm gpuinfo_*'
  call.call_simple_shell(work_dir, cmd_2)

  return device, usage
