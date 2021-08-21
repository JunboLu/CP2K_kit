#! /usr/env/bin python

from CP2K_kit.tools import data_op

def log_logo():

  logo = '''
********************************************************************************
********************************************************************************
**                                                                            **
**                                                                            **
**       ######    #####    ######  ##   ##       ##   ##  ##  ########       **
**     ##        ##     ##      ##  ##  ##        ##  ##   ##     ##          **
**     ##        ##     ##      ##  ## ##         ## ##    ##     ##          **
**     ##        #######    ######  ####   ###### ####     ##     ##          **
**     ##        ##         ##      ## ##         ## ##    ##     ##          **
**       ######  ##         ##      ##  ##        ##  ##   ##     ##          **
**               ##         ######  ##   ##       ##   ##  ##     ##          **
**                                                                            **
**                        CP2K-KIT : A helper of CP2K                         **
**                        Email : lujunbo15@gmail.com                         **
**                                                                            **
********************************************************************************
'''

  print (logo, flush=True)

def log_error(error, error_lable='Error'):

  '''
  log_error : print error information

  Args :
    error : string
      error is the error information
    error_lable : string
      error_lable is the lable of error
  Returns :
    none
  '''

  print (error_lable.center(80,'*'), flush=True)
  print (data_op.str_wrap(error,80), flush=True)
  print ('\n', flush=True)
