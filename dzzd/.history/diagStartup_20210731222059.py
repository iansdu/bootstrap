# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 18:55:01 2019

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
from scipy.optimize import fmin
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH,'c:\\Program Files (x86)\\REFPROP\\')
CP.get_global_param_string("REFPROP_version")
from FuncMod import *
from AF_Pro import AF
import csv
from DPMod import * #DPMod cal #initial parameters
from ODMod import *
import xlrd
from HPC_GTMap import HPC_GTMap
from HPT_GTMap import HPT_GTMap



def Diagiter(xone):

    DEDG.DGC =xone[0]
    DEDG.DEC =xone[1]
    DEDG.DGT =xone[2]
    DEDG.DET =xone[3]

    error = numpy.zeros(4)
    #try:
    SV25.P = SV25_P
    SV3.P = SV3_P
    SV25.T = SV25_T
    HCom_DPIP.PR_OD =SV3.P/SV25.P
    HPC_OD(SV25,SV3,HCom_DPIP,3,SV251,SV252)
    #Divide_OD(SV251,SV_HPTRotor,SV_LPTRotor,Rotor_Divide)
    #Divide_OD(SV252,SV_HPTNGV,SV_LPTNGV,NGV_Divide)
    SVF.MF = SVF_MF
    Burn_OD(SV3,SV4,SVF,Combustor_DPIP)
    #Mix_OD(SV4,SV41,SV_HPTNGV,HPTNGV_DP,SVF) #PD = 1,EF delete
    #SV41midMF = SV41.MF
    HPT_DPIP.PR_OD = SV4.P/(SV45_P/(1-PD_HPout_DPIP.PD_DP))
    HPT_OD(SV4,SV45,HPT_DPIP)
    #Mix_OD(SV42,SV43,SV_HPTRotor,HPTRotor_DP,SVF)
    #PD_OD(SV42,SV45,PD_HPout_DPIP)

    
    error[0] = abs((HCom_DPIP.CW-HPT_DPIP.CW)/HCom_DPIP.CW)
    error[1] = abs(((SV4.MF-SV45.MF)/SV4.MF))
    error[2] = abs((SV3_T-SV3.T)/SV3_T)
    error[3] = abs((SV45_T-SV45.T)/SV45_T)
    #except:
    #    error = numpy.array([100,100,100,100])
    #    pass
    print(error.sum())
    return(error)
# field data assignment

SV1 = AF()
SV2 = AF()
SV21 = AF()
SV24 = AF()
SV24B = AF()   
SV25 = AF()
SV251 = AF()
SV252 = AF()
SV3 = AF()
SV4 = AF()
SV41 = AF()
SV42 = AF()
SV43 = AF()
SV44 = AF()
SV45 = AF()
SV46 = AF()
SV47 = AF()
SV48 = AF()   
SV49 = AF()
SV5 = AF()
SV6 = AF()
SV7 = AF()
SV9 = AF()
SV11 = AF()
SV13 = AF()
SV14 = AF()
SV14B = AF()      
SV15 = AF()
SV16 = AF()
SV17 = AF()
SV19 = AF()    
SV_Ecs = AF()
SV_Hand = AF()    
SV_TBV = AF()
SV_VBV = AF()  
SV_Byp = AF()      
SV_HPTRotor = AF()
SV_HPTNGV = AF()
SV_HPTExit = AF()
SV_LPTRotor = AF()
SV_LPTNGV = AF()
SV_LPTExit = AF()
SVF = AF()


book = xlrd.open_workbook('input.xlsx')
sheet1 = book.sheets()[0]
nrows = sheet1.nrows



fielddata_NH = sheet1.col_values(0)
fielddata_MFF = sheet1.col_values(1)#kg/s
fielddata_T25 = sheet1.col_values(2)+273.15
fielddata_P25 = sheet1.col_values(3)
fielddata_T3 = sheet1.col_values(4) +273.15
fielddata_P3 = sheet1.col_values(5)
fielddata_T45 = sheet1.col_values(6)+273.15
fielddata_P45 = sheet1.col_values(7)
fielddata_HPCB1_valve = 1
fielddata_HPCB2_valve = 1
fielddata_HPTRotorDivide_valve = 1
fielddata_HPTNGVDivide_valve = 1
diag = numpy.zeros((nrows,4))
mapline = numpy.zeros((nrows,8))
for i in range(nrows):
    HCom_DPIP.N = fielddata_NH [i]
    HPT_DPIP.N = fielddata_NH [i]
    SV25_T = fielddata_T25 [i]+273.15
    SV25_P = fielddata_P25 [i]*1000
    SV3_T = fielddata_T3 [i]+273.15
    SV3_P = fielddata_P3 [i]*1000
    SV45_T = fielddata_T45 [i]+273.15
    SV45_P = fielddata_P45 [i]*1000
    SVF_MF = fielddata_MFF[i]
    HCom_DPIP.valve1 = 1
    HCom_DPIP.valve2 = 1
    Rotor_Divide.valve = 1
    NGV_Divide.valve = 1

    SFguess = [0,0,0,0]
    Finalresult = fsolve(Diagiter, SFguess)
    SFout = Finalresult
    diag[i] = SFout

    mapline[i] = [
                  HCom_DPIP.CN_OD,HCom_DPIP.PR_AC,HCom_DPIP.CM_AC,HCom_DPIP.EF_AC,
                  HPT_DPIP.CN_OD,HPT_DPIP.PR_AC,HPT_DPIP.CM_AC,HPT_DPIP.EF_AC]
    numpy.savetxt('diag.csv',diag,delimiter=',')
    #numpy.savetxt('mapline.csv',mapline,delimiter=',')
    
