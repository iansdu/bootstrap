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

def Comp_DP(MF_in,T_in, P_in,T_out,P_out):    
    global Comp_PR_DP,Comp_EFF_DP,Comp_CW_DP
    Comp_PR_DP = P_out/P_in
    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
    fluid.set_mass_fractions([0.2314, 0.7552, 0.0129,0.0005])
    frac_Comp = fluid.get_mole_fractions()
    
    H_in = PropsSI('H','T',T_in,'P',P_in,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    S_in = PropsSI('S','T',T_in,'P',P_in,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_oI = PropsSI('H','S',S_in,'P',P_out,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))   
    H_out = PropsSI('H','P',P_out,'T',T_out,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    Comp_EFF_DP = (H_oI-H_in)/(H_out-H_in)
    SW = H_out - H_in #comp specific work
    CW = MF_in * SW #comp total work
    #bleed one
    PR_B1 = Comp_PRB1 * Comp_PR_DP
    P_B1 = P_in * PR_B1 
    MF_B1 = Comp_WACB1 * MF_in  
    ETA_B1 = math.pow(Comp_EFF_DP,Comp_PRB1)
    H_1I = PropsSI('H','P',P_B1,'S',S_in,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_B1 = H_in + (H_1I - H_in) / ETA_B1
    T_B1 = PropsSI('T','P',P_B1,'H',H_B1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    CW_B1M = (H_out - H_B1) * MF_B1

    #bleed two
    PR_B2 = Comp_PRB2 * Comp_PR_DP
    P_B2 = P_in * PR_B2
    MF_B2 = Comp_WACB2 * (MF_in - MF_B1)
    ETA_B2 = math.pow(Comp_EFF_DP,Comp_PRB2)
    H_2I = PropsSI('H','P',P_B2,'S',S_in,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_B2 = H_in + (H_2I - H_in) / ETA_B2
    T_B2 = PropsSI('T','P',P_B2,'H',H_B2,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    CW_B2M = (H_out - H_B2) * MF_B2   
    #MFin
    MF_out = (1-Comp_WACB1-Comp_WACB2+Comp_WACB1*Comp_WACB2)*MF_in 
    Comp_CW_DP = CW - CW_B1M - CW_B2M
    return MF_out,H_out,MF_B1, P_B1, T_B1, H_B1, MF_B2, P_B2, T_B2, H_B2

