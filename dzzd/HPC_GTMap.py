# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:54:42 2019

@author: zhaoqi
"""

import numpy
import math
import scipy
import csv

datadir = 'map'

class HPC_GTMap():
    # with open(datadir+'/HPCMap.csv', 'r') as f:
    with open('C:/Users/ian_s/Desktop/Test/dzzd/map/HPCMap.csv', 'r') as f:
        reader = csv.reader(f)
        result = list(reader)
        result = numpy.array(result).astype(numpy.float64).tolist()
    nCN = int(len(result[0]))
    nPOINTS = int((len(result)-1)/3)
    AX = numpy.array(result[0]) #speed line, CN or PCN???

    #PR
    BX = numpy.mat(numpy.array(result[1:nPOINTS+1]).T.tolist())
    # MASS
    CX = numpy.mat(numpy.array(result[nPOINTS+1:nPOINTS*2+1]).T.tolist())

    #Eff
    DX = numpy.mat(numpy.array(result[nPOINTS*2+1:nPOINTS*3+1]).T.tolist())