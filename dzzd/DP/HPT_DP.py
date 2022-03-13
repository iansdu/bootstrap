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

def HPT_DP(MF_in, T_in, P_in, H_in, PR):
    global HPT_TW_DP,HPT_PR_DP, PRHT_SF_DP, CMHT_SF_DP, ETAHT_SF_DP, HPT_N_ref, HPT_EFF_DP
    HPT_N_ref = Comp_N_DP / math.sqrt(T_in / TSL)
    HPT_SQRTA = math.sqrt(T_in / TSL)
    CM_in = MF_in * HPT_SQRTA / (P_in / PSL)
    HPT_TW_DP = Comp_CW_DP
    P_ou = P_in/PR
    H_ou = H_in - HPT_TW_DP / MF_in
    T_ou = PropsSI('T','H',H_ou,'P',P_ou,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    S_in = PropsSI('S','T',T_in,'P',P_in,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    H_oI = PropsSI('H','P',P_ou,'S',S_in,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    HPT_EFF_DP = (H_in - H_ou)/(H_in - H_oI)
    (HPT_PR_R, HPT_CM_R, HPT_ETA_R) = ReadMap(HPT_CN_DP,HPT_PR_MAP, 0, LM2500Map.LM2500HPT)
    PRHT_SF_DP = (HPT_PR_DP - 1.0) / (HPT_PR_R - 1.0)
    CMHT_SF_DP = CM_in / HPT_CM_R
    ETAHT_SF_DP = HPT_EFF_DP / HPT_ETA_R
    MF_ou = MF_in  
    return MF_ou,T_ou,P_ou, H_ou