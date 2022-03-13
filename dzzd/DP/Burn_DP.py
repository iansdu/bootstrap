# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:29 2019

@author: dell
"""
from CoolProp.CoolProp import PropsSI
import numpy
from scipy import interpolate
import math
import CoolProp
import scipy
import CoolProp.CoolProp as CP
import time
import datetime
from scipy.optimize import fsolve
from scipy.optimize import newton
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH,'c:\\Program Files (x86)\\REFPROP\\')
CP.get_global_param_string("REFPROP_version")

def Burn_DP(MF_in, T_in, P_in, H_in):
    global Burn_DP_PD, FuelF, frac, LHV, Burn_Drop_DP, ETABSF_DP, Burn_DP_Eff
    global c1, c2, c3, c4, ic4, c5, ic5, c6, N2, co2

    P_ou = P_in * (1.0 - Burn_PDamage_DP)
    Burn_Drop_DP = Burn_PDamage_DP * P_in / (MF_in * MF_in * T_in)

    H_ou = (MF_in * H_in + FuelFlow_DP * LHV + 87037 * FuelFlow_DP) * Burn_EFF_DP/(MF_in + FuelFlow_DP) 

    HD = H_ou - H_in

    #FAR = FuelFlow_DP / MF_in #FAR
    GasProp(MF_in, FuelFlow_DP)
    
    T_ou = PropsSI('T','H',H_ou,'P',P_ou,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    HT = T_ou - T_in
    MF_out = MF_in +FuelFlow_DP
    (HT, Burn_CM_R, Burn_ETA_R) = ReadMap(0.00001 * P_in, HT, 0, LM2500Map.FixedBURNERMap1)

    ETABSF_DP = Burn_EFF_DP / Burn_ETA_R
    return MF_out,T_ou, P_ou, H_ou
