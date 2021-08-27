#! /usr/env/bin python

import copy
import linecache
from collections import OrderedDict
from CP2K_kit.tools import data_op
from CP2K_kit.tools import log_info

#This module is used to read input file. It is a bit complicated!
#We have to use it now. But the further revisions are on the way.
#Currently, we haven't written enough comments. We hope it could
#work for a long time.

def get_keyword(part_range, inp):

  keyword = []
  keyword_index = []
  for i in range(part_range[0],part_range[1]+1,1):
    line_i = linecache.getline(inp,i)
    if ('&' in line_i):
      index_1 = line_i.index('&')
      index_2 = line_i.index('\n')
      keyword_index.append(i)
      keyword.append(line_i[index_1+1:index_2])

  return keyword, keyword_index

def get_keyword_block(keyword, keyword_index):

  keyword_block = []
  keyword_block_index = []
  for i in range(len(keyword)):
    key_same_num = 0
    end_key_same_num = 0
    for j in range(i,len(keyword),1):
      if ( keyword[i] == keyword[j] ):
        key_same_num = key_same_num + 1
      if ( keyword[i] != keyword[j] and keyword[i] in keyword[j]):
        end_key_same_num = end_key_same_num + 1
        if ( key_same_num == 0 or key_same_num == end_key_same_num ):
          keyword_block.append(keyword[i])
          keyword_block_index.append([keyword_index[i]+1,keyword_index[j]-1])
          break

  return keyword_block, keyword_block_index

def get_dump(keyword_block, keyword_block_index, inp):

  dump_dic = OrderedDict()

  sys_num = 0
  for i in range(len(keyword_block)):
    if ( keyword_block[i] == 'group' ):
      keyword_block[i] = 'group' + str(sys_num)
      sys_num = sys_num + 1

  sys_num = 0
  for i in range(len(keyword_block)):
    if ( keyword_block[i] == 'system' ):
      keyword_block[i] = 'system' + str(sys_num)
      sys_num = sys_num + 1

  sys_num = 0
  for i in range(len(keyword_block)):
    if ( keyword_block[i] == 'connect' ):
      keyword_block[i] = 'connect' + str(sys_num)
      sys_num = sys_num + 1

  if ( len(keyword_block) == 1 ):
    for i in range(keyword_block_index[0][0],keyword_block_index[0][1]+1,1):
      line_i = linecache.getline(inp, i)
      line_i_split = data_op.str_split(line_i, ' ')
      if (line_i_split[len(line_i_split)-1] == '\n'):
        line_i_split.remove('\n')
      if ( len(line_i_split) == 2 ):
        dump_dic[line_i_split[0]] = line_i_split[1].strip('\n')
      elif ( len(line_i_split) > 2 ):
        line_i_split[len(line_i_split)-1] = line_i_split[len(line_i_split)-1].strip('\n')
        dump_dic[line_i_split[0]] = line_i_split[1:len(line_i_split)]
  else:
    for i in range(0, len(keyword_block), 1):

      i_list = data_op.gen_list(keyword_block_index[i][0],keyword_block_index[i][1],1)
      i_dic = OrderedDict()
      sub_key = []
      sub_key_index = []

      for j in range(1, len(keyword_block), 1):
        if ( i != j ):
          if ( keyword_block_index[i][0] < keyword_block_index[j][0] \
               and keyword_block_index[i][1] > keyword_block_index[j][1] ):
            j_list = data_op.gen_list(keyword_block_index[j][0],keyword_block_index[j][1],1)
            sub_key.append(keyword_block[j])
            sub_key_index.append(j_list)
            j_list_temp = data_op.gen_list(keyword_block_index[j][0]-1,keyword_block_index[j][1]+1,1)
            i_list = [x for x in i_list if x not in j_list_temp]
          elif ( keyword_block_index[i][0] > keyword_block_index[j][0] \
                 and keyword_block_index[i][1] < keyword_block_index[j][1] ):
            i_list = []
            break

      if ( i_list == [] and sub_key != [] and i != 0 ):
        for j in range(len(sub_key)):
          j_dic = OrderedDict()
          for k in sub_key_index[j]:
            line_k = linecache.getline(inp, k)
            line_k_split = data_op.str_split(line_k, ' ')
            if (line_k_split[len(line_k_split)-1] == '\n'):
              line_k_split.remove('\n')
            if ( len(line_k_split) == 2 ):
              j_dic[line_k_split[0]] = line_k_split[1].strip('\n')
            elif ( len(line_k_split) > 2 ):
              line_k_split[len(line_k_split)-1] = line_k_split[len(line_k_split)-1].strip('\n')
              j_dic[line_k_split[0]] = line_k_split[1:len(line_k_split)]
          i_dic[sub_key[j]] = j_dic
        dump_dic[keyword_block[i]] = i_dic

      elif ( i_list != [] and sub_key != [] and i != 0 ):
        for j in i_list:
          line_j = linecache.getline(inp, j)
          line_j_split = data_op.str_split(line_j, ' ')
          if (line_j_split[len(line_j_split)-1] == '\n'):
            line_j_split.remove('\n')
          if ( len(line_j_split) == 2 ):
            i_dic[line_j_split[0]] = line_j_split[1].strip('\n')
          elif ( len(line_j_split) > 2 ):
            line_j_split[len(line_j_split)-1] = line_j_split[len(line_j_split)-1].strip('\n')
            i_dic[line_j_split[0]] = line_j_split[1:len(line_j_split)]
        for j in range(len(sub_key)):
          if ( sub_key[j] in dump_dic ):
            dump_dic.pop(sub_key[j])
          j_dic = OrderedDict()
          for k in sub_key_index[j]:
            line_k = linecache.getline(inp, k)
            line_k_split = data_op.str_split(line_k, ' ')
            if (line_k_split[len(line_k_split)-1] == '\n'):
              line_k_split.remove('\n')
            if ( len(line_k_split) == 2 ):
              j_dic[line_k_split[0]] = line_k_split[1].strip('\n')
            elif ( len(line_k_split) > 2 ):
              line_k_split[len(line_k_split)-1] = line_k_split[len(line_k_split)-1].strip('\n')
              j_dic[line_k_split[0]] = line_k_split[1:len(line_k_split)]
          i_dic[sub_key[j]] = j_dic
        dump_dic[keyword_block[i]] = i_dic

      elif ( i_list != [] and sub_key == []):
        for j in i_list:
          line_j = linecache.getline(inp, j)
          line_j_split = data_op.str_split(line_j, ' ')
          if (line_j_split[len(line_j_split)-1] == '\n'):
            line_j_split.remove('\n')
          if ( len(line_j_split) == 2 ):
            i_dic[line_j_split[0]] = line_j_split[1].strip('\n')
          elif ( len(line_j_split) > 2 ):
            line_j_split[len(line_j_split)-1] = line_j_split[len(line_j_split)-1].strip('\n')
            i_dic[line_j_split[0]] = line_j_split[1:len(line_j_split)]
        dump_dic[keyword_block[i]] = i_dic

      elif ( i_list != [] and sub_key != [] and i == 0):
        for j in i_list:
          line_j = linecache.getline(inp, j)
          line_j_split = data_op.str_split(line_j, ' ')
          if (line_j_split[len(line_j_split)-1] == '\n'):
            line_j_split.remove('\n')
          if ( len(line_j_split) == 2 ):
            dump_dic[line_j_split[0]] = line_j_split[1].strip('\n')
          elif ( len(line_j_split) > 2 ):
            line_j_split[len(line_j_split)-1] = line_j_split[len(line_j_split)-1].strip('\n')
            dump_dic[line_j_split[0]] = line_j_split[1:len(line_j_split)]
        for j in range(len(sub_key)):
          j_dic = OrderedDict()
          for k in sub_key_index[j]:
            line_k = linecache.getline(inp, k)
            line_k_split = data_op.str_split(line_k, ' ')
            if (line_k_split[len(line_k_split)-1] == '\n'):
              line_k_split.remove('\n')
            if ( len(line_k_split) == 2 ):
              j_dic[line_k_split[0]] = line_k_split[1].strip('\n')
            elif ( len(line_k_split) > 2 ):
              line_k_split[len(line_k_split)-1] = line_k_split[len(line_k_split)-1].strip('\n')
              j_dic[line_k_split[0]] = line_k_split[1:len(line_k_split)]
          dump_dic[sub_key[j]] = j_dic

  return dump_dic

def dump_info(work_dir, inp_file, f_key):

  f_key_copy = copy.deepcopy(f_key)
  input_file = work_dir + '/' + inp_file
  whole_line_num = len(open(input_file).readlines())
  inp_parse = []

  f_key_range = []

  for keyword in f_key_copy:
    keyword_range = []
    for i in range(whole_line_num):
      line_i = linecache.getline(input_file, i+1)
      if ( keyword in line_i and '&' in line_i ):
        keyword_range.append(i+1)
    f_key_range.append(keyword_range)

  for i in range(len(f_key_range)):
    if ( len(f_key_range[i]) != 2):
      log_info.log_error('The %s parse is incompleted' %(f_key[i]))
      exit()

  for i in range(len(f_key_copy)):
    f_key_i, f_key_i_index = get_keyword(f_key_range[i], input_file)

    f_key_i_block, f_key_i_block_index = \
    get_keyword_block(f_key_i, f_key_i_index)

    f_key_i_dic = get_dump(f_key_i_block, f_key_i_block_index, input_file)

    inp_parse.append(f_key_i_dic)

  inp_parse_copy = copy.deepcopy(inp_parse)

  for i in range(len(inp_parse)):
    for key1 in inp_parse[i]:
      for key2 in inp_parse[i]:
        if ( type(inp_parse[i][key2]) == OrderedDict and key1 in inp_parse[i][key2].keys() ):
          inp_parse_copy[i].pop(key1)

  return inp_parse_copy

if __name__ == '__main__':
  from CP2K_kit.analyze import read_input
  work_dir = '/home/lujunbo/code/github/CP2K_kit/example/analyze/coord_num'
  inp_file = 'input.inp'
  f_key = ['geometry']
  inp_parse = read_input.dump_info(work_dir, inp_file, f_key)
  print (inp_parse)
