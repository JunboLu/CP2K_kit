#! /usr/env/bin python

import os
import math
from CP2K_kit.deepff import process
from CP2K_kit.tools import *

def run_cp2kfrc(work_dir, iter_id, cp2k_exe, parallel_exe, cp2k_env_file, cp2k_job_per_node, proc_num_per_node, host, ssh, atoms_num_tot):

  '''
  run_force: perform cp2k force calculation

  Args:
    work_dir: string
      work_dir is the workding directory of CP2K_kit.
    iter_id: int
      iter_id is current iteration number.
    cp2k_exe: string
      cp2k_exe is the cp2k executable file.
    parallel_exe: string
      parallel_exe is the parallel executable file.
    cp2k_env_file: string
      cp2k_env_file is the environment setting file of cp2k.
    cp2k_job_per_node: int
      cp2k_job_per_node is the job number of cp2k in each node.
    proc_num_per_node: 1-d int list
      proc_num_per_node is the numbers of processor in each node.
    host: 1-d string list
      host is the name of computational nodes.
    ssh: bool
      ssh is whether we need to ssh.
    atoms_num_tot: 1-d dictionary, dim = num of lammps systems
      example: {1:192,2:90}
  Returns:
    none
  '''

  import subprocess

  cp2k_calc_dir = ''.join((work_dir, '/iter_', str(iter_id), '/03.cp2k_calc'))
  sys_num = process.get_sys_num(cp2k_calc_dir)

  #check generating cp2k tasks
  check_cp2k_gen = []
  for i in range(sys_num):
   cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
   task_num, task_dir = process.get_task_num(cp2k_sys_dir, True)
   for j in range(task_num):
     cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/', task_dir[j]))
     traj_num = process.get_traj_num(cp2k_sys_task_dir)
     for k in range(traj_num):
       cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(k)))
       inp_file_name_abs = ''.join((cp2k_sys_task_traj_dir, '/input.inp'))
       if ( os.path.exists(inp_file_name_abs) and os.path.getsize(inp_file_name_abs) != 0 ):
         check_cp2k_gen.append(0)
       else:
         check_cp2k_gen.append(1)
  if ( all(i == 0 for i in check_cp2k_gen) ):
    str_print = 'Success: generate cp2k tasks in %s' %(cp2k_calc_dir)
    str_print = data_op.str_wrap(str_print, 80, '  ')
    print (str_print, flush=True)
  else:
    log_info.log_error('Generating cp2k tasks error, please check iteration %d' %(iter_id))
    exit()

  #run cp2k tasks
  for i in range(sys_num):
    cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
    task_num, task_dir = process.get_task_num(cp2k_sys_dir, True)
    host_name_proc = []
    for l in range(len(host)):
      host_name_proc.append(''.join((str(proc_num_per_node[l]), '/', host[l])))
    host_info = data_op.comb_list_2_str(host_name_proc, ',')

    for j in range(task_num):
      cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/', task_dir[j]))
      traj_num = process.get_traj_num(cp2k_sys_task_dir)
      calculated_id = 0
      for k in range(traj_num):
        cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(k)))

        frc_file=''.join((cp2k_sys_task_traj_dir, '/cp2k-1_0.xyz'))
        if ( os.path.exists(frc_file) and len(open(frc_file, 'r').readlines()) == atoms_num_tot[i]+5 ):
          calculated_id = k
        else:
          break

      if ( calculated_id != traj_num-1 ):
        run_start = calculated_id
        run_end = run_start+cp2k_job_per_node*len(host)-1
        if ( run_end > traj_num-1 ):
          run_end=traj_num-1
        cycle = math.ceil((traj_num-run_start)/(cp2k_job_per_node*len(host)))

        for k in range(cycle):
          tot_mpi_num_list = []
          for proc_num in proc_num_per_node:
            mpi_num_list = data_op.int_split(proc_num, cp2k_job_per_node)
            for l in range(len(mpi_num_list)):
              if ( mpi_num_list[l]%2 != 0 and mpi_num_list[l]>1 ):
                mpi_num_list[l] = mpi_num_list[l]-1
            tot_mpi_num_list.append(mpi_num_list)
          tot_mpi_num_list = data_op.list_reshape(tot_mpi_num_list)[0:(run_end-run_start+1)]
          mpi_num_str = data_op.comb_list_2_str(tot_mpi_num_list, ' ')
          task_job_list = data_op.gen_list(run_start, run_end, 1)
          task_job_str = data_op.comb_list_2_str(task_job_list, ' ')

          run_1 = '''
#! /bin/bash
task_job="%s"
mpi_num="%s"
direc=%s
parallel_exe=%s
task_job_arr=(${task_job///})
mpi_num_arr=(${mpi_num///})
num=${#task_job_arr[*]}
for ((i=0;i<=num-1;i++));
do
task_job_mpi_num_arr[i]="${task_job_arr[i]} ${mpi_num_arr[i]}"
done
''' %(task_job_str, mpi_num_str, cp2k_sys_task_dir, parallel_exe)
          if ssh:
            run_2 = '''
for i in "${task_job_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d -S %s --controlmaster --sshdelay 0.2 $direc/produce.sh {} $direc
''' %( cp2k_job_per_node, host_info)
          else:
            run_2 = '''
for i in "${task_job_mpi_num_arr[@]}"; do echo "$i"; done | $parallel_exe -j %d --delay 0.2 $direc/produce.sh {} $direc
''' %( cp2k_job_per_node)

          produce = '''
#! /bin/bash
source %s
x=$1
direc=$2
x_arr=(${x///})
new_direc=$direc/traj_${x_arr[0]}
cd $new_direc
if [ -f "cp2k-1_0.xyz" ]; then
rm cp2k-1_0.xyz
fi
mpirun -np ${x_arr[1]} %s $new_direc/input.inp 1> $new_direc/cp2k.out 2> $new_direc/cp2k.err
converge_info=`grep "SCF run NOT converged" cp2k.out`
if [ $? -eq 0 ]; then
wfn_line=`grep -n "WFN_RESTART_FILE_NAME" input.inp`
if [ $? -eq 0 ]; then
line=`grep -n "WFN_RESTART_FILE_NAME" input.inp | awk -F ":" '{print $1}'`
sed -i ''$line's/.*/    WFN_RESTART_FILE_NAME .\/cp2k-RESTART.wfn/' input.inp
else
line=`grep -n "POTENTIAL_FILE_NAME" input.inp | awk -F ":" '{print $1}'`
sed -i ''$line' s/^/    WFN_RESTART_FILE_NAME .\/cp2k-RESTART.wfn\\n/' input.inp
fi
if [ -f "cp2k-1_0.xyz" ]; then
rm cp2k-1_0.xyz
fi
mpirun -np ${x_arr[1]} %s $new_direc/input.inp 1> $new_direc/cp2k.out 2> $new_direc/cp2k.err
fi
cd %s
''' %(cp2k_env_file, cp2k_exe, cp2k_exe, work_dir)

          run_file_name_abs = ''.join((cp2k_sys_task_dir, '/run.sh'))
          with open(run_file_name_abs, 'w') as f:
            f.write(run_1+run_2)

          produce_file_name_abs = ''.join((cp2k_sys_task_dir, '/produce.sh'))
          with open(produce_file_name_abs, 'w') as f:
            f.write(produce)

          subprocess.run('chmod +x run.sh', cwd=cp2k_sys_task_dir, shell=True)
          subprocess.run('chmod +x produce.sh', cwd=cp2k_sys_task_dir, shell=True)
          try:
            subprocess.run("bash -c './run.sh'", cwd=cp2k_sys_task_dir, shell=True)
          except subprocess.CalledProcessError as err:
            log_info.log_error('Running error: %s command running error in %s' %(err.cmd, cp2k_sys_task_dir))

          run_start = run_start + cp2k_job_per_node*len(host)
          run_end = run_end + cp2k_job_per_node*len(host)
          if ( run_end > traj_num-1):
            run_end = traj_num-1
      else:
        cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(traj_num-1)))
        frc_file_name_abs = ''.join((cp2k_sys_task_traj_dir, '/cp2k-1_0.xyz'))
        if ( os.path.exists(frc_file_name_abs) ):
          cmd = "rm %s" %(frc_file_name_abs)
          call.call_simple_shell(cp2k_sys_task_traj_dir, cmd)
        cmd = "bash -c 'source %s && mpirun -np %d %s input.inp 1> cp2k.out 2> cp2k.err'" \
              %(cp2k_env_file, sum(proc_num_per_node), cp2k_exe)
        try:
          subprocess.run(cmd, cwd=cp2k_sys_task_traj_dir, shell=True)
        except subprocess.CalledProcessError as err:
          log_info.log_error('Running error: %s command running error in %s' %(err.cmd, cp2k_sys_task_traj_dir))

  #check running cp2k tasks
  check_cp2k_run = []
  for i in range(sys_num):
    cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
    task_num = process.get_task_num(cp2k_sys_dir)
    check_cp2k_run_i = []
    for j in range(task_num):
      check_cp2k_run_ij = []
      cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(j)))
      traj_num = process.get_traj_num(cp2k_sys_task_dir)
      for k in range(traj_num):
        cp2k_sys_task_traj_dir = ''.join((cp2k_sys_task_dir, '/traj_', str(k)))
        frc_file_name_abs = ''.join((cp2k_sys_task_traj_dir, '/cp2k-1_0.xyz'))
        if ( os.path.exists(frc_file_name_abs) and len(open(frc_file_name_abs, 'r').readlines()) == atoms_num_tot[i]+5 ):
          check_cp2k_run_ij.append(0)
        else:
          check_cp2k_run_ij.append(1)
      check_cp2k_run_i.append(check_cp2k_run_ij)
    check_cp2k_run.append(check_cp2k_run_i)
  if ( all(i == 0 for i in data_op.list_reshape(data_op.list_reshape(check_cp2k_run))) ):
    print ('  Success: ab initio force calculations for %d systems by cp2k' %(sys_num), flush=True)
  else:
    for i in range(sys_num):
      cp2k_sys_dir = ''.join((cp2k_calc_dir, '/sys_', str(i)))
      task_num = process.get_task_num(cp2k_sys_dir)
      for j in range(task_num):
        cp2k_sys_task_dir = ''.join((cp2k_sys_dir, '/task_', str(j)))
        traj_num = process.get_traj_num(cp2k_sys_task_dir)
        failure_task_id = [index for (index,value) in enumerate(check_cp2k_run[i][j]) if value==1]
        if ( len(failure_task_id) != traj_num and len(failure_task_id) != 0 ):
          failure_task_id_str = data_op.comb_list_2_str(failure_task_id, ' ')
          str_print = '  Warning: ab initio force calculations fail for tasks %s in system %d in task %d by cp2k' %(failure_task_id_str, i, j)
          str_print = data_op.str_wrap(str_print, 80, '  ')
          print (str_print, flush=True)
        elif ( len(failure_task_id) == traj_num):
          log_info.log_error('Running error: ab initio force calculations running error, please check iteration %d' %(iter_id))
          exit()

if __name__ == '__main__':
  from collections import OrderedDict
  from CP2K_kit.tools import read_input
  from CP2K_kit.deepff import cp2k_run
  from CP2K_kit.deepff import check_deepff

  work_dir = '/home/lujunbo/WORK/Deepmd/CP2K_kit/co2/md_mtd'
  deepff_key = ['deepmd', 'lammps', 'cp2k', 'model_devi', 'environ']
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  read_input.dump_info(work_dir, 'input.inp', deepff_key)
  proc_num = 4
  deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic = \
  check_deepff.check_inp(deepmd_dic, lammps_dic, cp2k_dic, model_devi_dic, environ_dic, proc_num)

  #Test gen_cp2kfrc_file function
  coord = [[['O','-9.64','-0.71','5.80'],['H','-10.39','-1.31','6.15'],['H','-8.89','-35.4','6.37']],
            [['O','-2.64','-7.14','5.52'],['H','-2.89','-6.23','5.10'],['H','-1.70','-7.36','5.28']]]
  box = [[['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']],
         [['10.0','0.0','0.0'],['0.0','10.0','0.0'],['0.0','0.0','10.0']]]
#  cp2k_run.gen_cp2kfrc_file(cp2k_dic, work_dir, 1, 0, coord, box)

  #Test gen_cp2k_task function
  atoms_type_multi_sys = {0: {'C': 1, 'O': 2}, 1: {'C': 1, 'O': 2}}
  atoms_num_tot = {0:3, 1:3}
  struct_index = OrderedDict([(0, OrderedDict([(0, [237, 264, 275, 291, 331, 367, 422])])), (1, OrderedDict([(0, [])]))])
  conv_new_data_num = 5
  choose_new_data_num_limit = 100
  cp2k_run.gen_cp2k_task(cp2k_dic, work_dir, 17, atoms_type_multi_sys, atoms_num_tot, \
                         struct_index, conv_new_data_num, choose_new_data_num_limit, False)

  #Test run_cp2kfrc function
#  cp2k_run.run_cp2kfrc(work_dir, 0, environ_dic)
