# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:29 2019

@author: dell
"""
import numpy
from scipy import interpolate
import math
import scipy

def Mix(MF_Main, T_Main, P_Main, H_Main, MF_Ble, T_Ble, P_Ble, H_Ble):
    global frac
    MF_ou = MF_Main + MF_Ble
    H_ou = (H_Main * MF_Main + H_Ble * MF_Ble)/MF_ou
    P_ou = P_Main
    GasProp(MF_Main - FuelFlow_DP + MF_Ble, FuelFlow_DP)
    T_ou = PropsSI('T','H',H_ou,'P',P_ou,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac)) #frac should be updated
    return MF_ou, P_ou, T_ou, H_ou