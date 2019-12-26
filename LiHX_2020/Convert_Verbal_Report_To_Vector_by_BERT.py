#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 12:45:30 2018

@author: lihuixian
"""


#bert-serving-start -model_dir /model/chinese_L-12_H-768_A-12/ -num_worker=4 -port 5555 -port_out 5556


from bert_serving.client import BertClient

bc = BertClient(ip='192.168.13.19',port=5555,port_out=5556)
import os
import numpy
import glob
readPath = '/Users/lihuixian/Documents/2018analysis/bert3/first10.2'
savePath = '/Users/lihuixian/Documents/2018analysis/bert3/vector_first10.2'

files=glob.glob('%s/*.txt' %readPath)
for path in files:
    filename=os.path.basename(path)[:-4]
    print(filename)
    f = open(path) 
    lines = f.readlines()  
    avector = bc.encode(lines)
    print(avector)
    numpy.savetxt(r'%s/%s.csv' %(savePath,filename),avector)
    print('写入成功')        