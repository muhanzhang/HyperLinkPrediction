import numpy as np
import scipy.io as sio
dataname = 'chuancai.txt'
with open('data/'+dataname) as f:
  for line in f:
    tmp = line.split(':')
    dish = tmp[0]
    food = tmp[1].split(',')

