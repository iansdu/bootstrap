# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:29 2019

@author: dell
"""
import numpy
from scipy import interpolate
import math
import scipy

def Loc(a_M,Val, NumT):
    Coeff = 0
    locV = 0
    for i in range(NumT):
        if a_M[i] > Val:
            locV = i;
            Coeff = (Val - a_M[i - 1]) / (a_M[i] - a_M[i - 1]);
            break;
    return Coeff, locV
