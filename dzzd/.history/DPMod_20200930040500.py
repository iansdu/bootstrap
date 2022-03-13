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
from DP_IP import *
from FuncMod import *
from AF_Pro import AF
import csv
from HPC_GTMap import HPC_GTMap
from HPT_GTMap import HPT_GTMap
from LPC_GTMap import LPC_GTMap
from LPT_GTMap import LPT_GTMap
from Fan_GTMap import Fan_GTMap
import xlrd
import pandas as pd

def SFDP():
    (HCom_DPIP.PR_MAP, HCom_DPIP.CM_MAP, HCom_DPIP.EF_MAP) = ScaleMap(CN_Map,HPC_GTMap)
    (HPT_DPIP.PR_MAP, HPT_DPIP.CM_MAP, HPT_DPIP.EF_MAP) = ScaleMap(CN_Map,HPT_GTMap)

    HCom_DPIP.PR_SFDP = (HCom_DPIP.PR_DP - 1.0) / (HCom_DPIP.PR_MAP - 1.0)
    HCom_DPIP.CM_SFDP = HCom_DPIP.CM_DP / HCom_DPIP.CM_MAP #cm = mf*qrt T /P 
    HCom_DPIP.EF_SFDP = HCom_DPIP.EF_DP / HCom_DPIP.EF_MAP
    
    HPT_DPIP.PR_SFDP = (HPT_DPIP.PR_DP - 1.0) / (HPT_DPIP.PR_MAP - 1.0)
    HPT_DPIP.CM_SFDP = HPT_DPIP.CM_DP / HPT_DPIP.CM_MAP #cm = mf*qrt T /P 
    HPT_DPIP.EF_SFDP = HPT_DPIP.EF_DP / HPT_DPIP.EF_MAP 
    zzz=1

Tamb = 273.15
Pamb = 101325
baseline = xlrd.open_workbook('baseline.xlsx')
baselinesheet = baseline.sheets()[0]
baseline_NL = baselinesheet.col_values(0)
baseline_NH = baselinesheet.col_values(1)
NH_DP  = baseline_NH[0]
NL_DP = baseline_NL[0]
CN_Map = 1

SV1_DP = AF()
SV2_DP = AF()
SV21_DP = AF()
SV24_DP = AF()
SV24B_DP = AF()   
SV25_DP = AF()
SV251_DP = AF()
SV252_DP = AF()
SV3_DP = AF()
SV4_DP = AF()
SV41_DP = AF()
SV42_DP = AF()
SV43_DP = AF()
SV44_DP = AF()
SV45_DP = AF()
SV46_DP = AF()
SV47_DP = AF()
SV48_DP = AF()   
SV49_DP = AF()
SV5_DP = AF()
SV6_DP = AF()
SV7_DP = AF()
SV9_DP = AF()
SV11_DP = AF()
SV13_DP = AF()
SV14_DP = AF()
SV14B_DP = AF()      
SV15_DP = AF()
SV16_DP = AF()
SV17_DP = AF()
SV19_DP = AF()    
SV_Ecs_DP = AF()
SV_Hand_DP = AF()    
SV_TBV_DP = AF()
SV_VBV_DP = AF()  
SV_Byp_DP = AF()      
SV_HPTRotor_DP = AF()
SV_HPTNGV_DP = AF()
SV_HPTExit_DP = AF()
SV_LPTRotor_DP = AF()
SV_LPTNGV_DP = AF()
SV_LPTExit_DP = AF()
SVF_DP = AF()

FanBypass_DPIP = FanBypass()
Mix_VBV_DPIP = Mix_DPIP()
PDtotal_ExCul_DPIP= PressureDamage()
Lost_ExCul_DP = Lost_DPIP()
PD_InCul_DPIP= PressureDamage()
Lost_ExLC_DP = Lost_DPIP()
Divide_VBV_DPIP = Lost_DPIP()
Fan_DPIP =  Com_DPIP()
LCom_DPIP = Com_DPIP()
HPTNGV_Lost_DP = Lost_DPIP()
HPTRotor_DP = Mix_DPIP()
HPTNGV_DP = Mix_DPIP()
HCom_DPIP = Com_DPIP()
Combustor_DPIP = Burn_DPIP()
HPT_DPIP = Turb_DPIP()
LPT_DPIP = Turb_DPIP()
PD_HPout_DPIP=PressureDamage()
PD_LPout_DPIP=PressureDamage()
LPTRotor_DP = Mix_DPIP()
LPTNGV_DP = Mix_DPIP()
Rotor_Divide = Lost_DPIP()
NGV_Divide = Lost_DPIP()
DEDG = DEDG()

DPIP_para = pd.read_csv('DPIP_para.csv',encoding='utf-8',header=None)#DPIP_para.values[53][0]

HCom_DPIP.CN_DP=DPIP_para.values[0][0]
HCom_DPIP.PR_DP=DPIP_para.values[1][0]
HCom_DPIP.CM_DP=DPIP_para.values[2][0]
HCom_DPIP.EF_DP=DPIP_para.values[3][0]
HCom_DPIP.N_ref=DPIP_para.values[4][0]
HCom_DPIP.PRB1=DPIP_para.values[5][0]
HCom_DPIP.WACB1=DPIP_para.values[6][0]
HCom_DPIP.EFB1=DPIP_para.values[7][0]
HCom_DPIP.CMCB1=DPIP_para.values[8][0]
HCom_DPIP.PRRB1=DPIP_para.values[9][0]
HCom_DPIP.PRB2=DPIP_para.values[10][0]
HCom_DPIP.WACB2=DPIP_para.values[11][0]
HCom_DPIP.EFB2=DPIP_para.values[12][0]
HCom_DPIP.CMCB2=DPIP_para.values[13][0]
HCom_DPIP.PRRB2=DPIP_para.values[14][0]
Rotor_Divide.WACB=DPIP_para.values[15][0]
Rotor_Divide.valve_DP=DPIP_para.values[16][0]
NGV_Divide.WACB=DPIP_para.values[17][0]
NGV_Divide.valve_DP=DPIP_para.values[18][0]
Combustor_DPIP.EF=DPIP_para.values[19][0]
Combustor_DPIP.FuelEnthalpy=DPIP_para.values[20][0]
Combustor_DPIP.LHV=DPIP_para.values[21][0]
Combustor_DPIP.PD_DP=DPIP_para.values[22][0]
Combustor_DPIP.PD_Ref=DPIP_para.values[23][0]
HPTNGV_DP.EF_DP=DPIP_para.values[24][0]
HPTNGV_DP.PD_DP=DPIP_para.values[25][0]
HPT_DPIP.CN_DP=DPIP_para.values[26][0]
HPT_DPIP.PR_DP=DPIP_para.values[27][0]
HPT_DPIP.CM_DP=DPIP_para.values[28][0]
HPT_DPIP.EF_DP=DPIP_para.values[29][0]
HPT_DPIP.N_ref=DPIP_para.values[30][0]
HPTRotor_DP.EF_DP=DPIP_para.values[31][0]
HPTRotor_DP.PD_DP=DPIP_para.values[32][0]
PD_HPout_DPIP.PD_DP=DPIP_para.values[33][0]

#center point calibration of maps
SFDP()











