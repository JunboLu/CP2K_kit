#! /usr/env/bin python

from CP2K_kit.tools import list_dic_op

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

def log_error(error):

  '''
  log_error : print error information

  Args :
    error : string
      error is the error information
  Returns :
    none
  '''

  print ('Error'.center(80,'*'), flush=True)
  print (list_dic_op.str_wrap(error,80), flush=True)

