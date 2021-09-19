#! /usr/env/bin python

import os
import platform
import subprocess
import linecache
import multiprocessing
from CP2K_kit.tools import call
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info

def get_host(work_dir):

  '''
  get_host: get hosts of nodes

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
   proc_num: int
     proc_num is the number of processors
   host: 1-d string list
     host is the name of node host
   ssh: bool
     ssh is whether we need to ssh nodes.
  '''

  env_dist = os.environ

  if 'PBS_NODEFILE' in env_dist or 'LSB_DJOB_HOSTFILE' in env_dist:
    if 'PBS_NODEFILE' in env_dist:
      node_info_file = env_dist['PBS_NODEFILE']
    elif 'LSB_DJOB_HOSTFILE' in env_dist:
      node_info_file = env_dist['LSB_DJOB_HOSTFILE']
    #Get proc_num
    cmd = "cat %s | wc -l" %(node_info_file)
    proc_num = int(call.call_returns_shell(work_dir, cmd)[0])
    #Get host name
    host = []
    proc_num_per_node = []
    cmd = "cat %s | sort | uniq -c" %(node_info_file)
    node_info = call.call_returns_shell(work_dir, cmd)
    for node in node_info:
      node_split = data_op.split_str(node, ' ')
      proc_num_per_node.append(int(node_split[0]))
      host.append(node_split[1])
    #Get ssh
    ssh = True

  elif 'SLURM_JOB_NODELIST' in env_dist:
    node_list = env_dist['SLURM_JOB_NODELIST']
    cmd = "scontrol show hostnames $SLURM_JOB_NODELIST"
    host = call.call_returns_shell(work_dir, cmd)
    cmd = "scontrol show node $SLURM_JOB_NODELIST | grep CPUAlloc"
    proc_num_per_node = []
    host_info = call.call_returns_shell(work_dir, cmd)
    for host_info_i in host_info:
      host_info_i_split = host_info_i.split('=')
      proc_num_per_node.append(int(host_info_i_split[1].split(' ')[0]))
    proc_num = int(env_dist['SLURM_NPROCS'])
    ssh = True
  else:
    proc_num = int(multiprocessing.cpu_count()/2)
    proc_num_per_node = [proc_num]
    host = [platform.node()]
    ssh = False

  return proc_num, proc_num_per_node, host, ssh

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
      line_1_split = data_op.split_str(line_1, ' ')
      device.append(line_1_split[1])
      line_2 = linecache.getline(gpuinfo_file, a_int+i*4+2)
      line_2_split = data_op.split_str(line_2, ' ')
      mem_used = float(line_2_split[8][0:line_2_split[8].index('M')])
      mem_tot = float(line_2_split[10][0:line_2_split[10].index('M')])
      usage.append(mem_used/mem_tot)

    linecache.clearcache()

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
    cmd = "touch %s" %(gpuinfo_file)
    call.call_simple_shell(work_dir, cmd)
    if ssh:
      check_gpu = '''
#! /bin/bash

ssh -T %s
nvidia-smi 1> /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
nvidia-smi >> %s
fi
exit
''' %(host[i], gpuinfo_file)

    else:
      check_gpu = '''
#! /bin/bash

nvidia-smi 1> /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
nvidia-smi >> %s
fi
''' %(gpuinfo_file)

    check_gpu_file = ''.join((work_dir, '/check_gpu.sh'))
    with open(check_gpu_file, 'w') as f:
      f.write(check_gpu)

    subprocess.run('chmod +x check_gpu.sh', cwd=work_dir, shell=True)
    try:
      subprocess.run("bash -c './check_gpu.sh'", cwd=work_dir, shell=True)
    except subprocess.CalledProcessError as err:
      log_info.log_error('Running error: %s command running error in %s' %(err.cmd, work_dir))

    device_i, usage_i = read_gpuinfo(work_dir, gpuinfo_file)
    device.append(device_i)
    usage.append(usage_i)

  cmd_1 = 'rm check_gpu.sh'
  call.call_simple_shell(work_dir, cmd_1)
  cmd_2 = 'rm gpuinfo_*'
  call.call_simple_shell(work_dir, cmd_2)

  return device, usage

def get_lmp_path(work_dir):

  '''
  get_lmp_path: get the path of lammps.

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    lmp_path: string
      lmp_path is the path of lammps.
  '''

  lmp_exe = call.call_returns_shell(work_dir, 'which lmp')
  if ( len(lmp_exe) == 0 or 'no lmp in' in lmp_exe[0] ):
    lmp_exe = call.call_returns_shell(work_dir, 'which lmp_mpi')
    if ( len(lmp_exe) == 0 or 'no lmp_mpi in' in lmp_exe[0] ):
      log_info.log_error('Envrionment error: can not find lmp executable file, please set the environment for lammps')
      exit()
  lmp_exe_split = data_op.split_str(lmp_exe[0], '/')
  lmp_path = data_op.comb_list_2_str(lmp_exe_split[:-2], '/', True)

  return lmp_exe[0], lmp_path

def get_mpi_path(work_dir):

  '''
  get_mpi_path: get the path of mpi

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    mpi_path: string
      mpi_path is the path of mpi.
  '''

  mpi_exe = call.call_returns_shell(work_dir, 'which mpirun')
  if ( len(mpi_exe) == 0 or 'no mpirun in' in mpi_exe[0] ):
    log_info.log_error('Envrionment error: can not find mpi executable file, please set the environment for mpi')
    exit()
  else:
    mpi_exe_split = data_op.split_str(mpi_exe[0], '/')
    mpi_path = data_op.comb_list_2_str(mpi_exe_split[:-2], '/', True)

  return lmp_exe, mpi_path

def get_dp_path(work_dir):

  '''
  get_dp_path: get the path of deepmd-kit

  Args:
    work_dir: string
      work_dir is the working directory of CP2K_kit.
  Returns:
    dp_path: string
      dp_path is the path of deepmd-kit.
  '''

  dp_exe = call.call_returns_shell(work_dir, 'which dp')
  if ( len(dp_exe) == 0 or 'no dp in' in dp_exe[0] ):
    log_info.log_error('Envrionment error: can not find dp executable file, please set the environment for deepmd-kit')
    exit()
  else:
    dp_exe_split = data_op.split_str(dp_exe[0], '/')
    dp_path = data_op.comb_list_2_str(dp_exe_split[:-2], '/', True)

  return dp_path
