#! /usr/env/bin python

import psutil
import subprocess
from CP2K_kit.tools import log_info
from subprocess import check_call, Popen, STDOUT, CalledProcessError

def call_simple_shell(exe_dir, cmd):

  '''
  call_simple_shell: call shell scripts, and we do not need returns.

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
      Example: '/home/lujunbo/code/github/CP2K_kit/tools'
    cmd: string
      cmd is the shell command.
      Example: 'ls -la'
  Returns:
    none
  '''

  #check_call is better. We could use call, check_call and run here.
  #But we need higher version python if we use run.
  try:
    check_call(cmd, cwd=exe_dir, shell=True)
  except CalledProcessError as err:
    log_info.log_error('Running error: %s command running error in %s' %(err.cmd, exe_dir))

def call_returns_shell(exe_dir, cmd):

  '''
  call_returns_shell: call shell scripts, and we need returns.

  Args:
    exe_dir: string
      exe_dir is the directory where shell script will be excuted.
      Example: '/home/lujunbo/code/github/CP2K_kit/tools'
    cmd: string
      cmd is the shell command.
      Example: 'ls -la'
  Returns:
    output: 1-d string list
  '''

  p = Popen(cmd, cwd=exe_dir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  (stdoutput, erroutput) = p.communicate()
  stdoutput = stdoutput.decode()
  output = stdoutput.split('\n')
  output.remove('')

  return output

def kills(parent_pid):

  '''
  kills: kill all process including parent and child processes.

  Args:
    parent_id: int
      parent_id is parent process.
  Returns:
    none
  '''

  parent = psutil.Process(parent_pid)
  for child in parent.children(recursive=True):
    child.kill()
  parent.kill()

if __name__ == '__main__':
  from CP2K_kit.tools import call
  #Test call_returns_shell function
  exe_dir = '/home/lujunbo/code/github/CP2K_kit/tools'
  cmd_1 = "grep -n 'savitzky_golay' numeric.py"
  out = call.call_returns_shell(exe_dir, cmd_1)
  print (out)

  #Test call_simple_shell function
  cmd_2 = 'mkdir kk'
  call.call_simple_shell(exe_dir, cmd_2)

