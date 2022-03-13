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
import xlrd
from HPC_GTMap import HPC_GTMap
from HPT_GTMap import HPT_GTMap



TSL = 273.15
PSL = 101325

#EE airinlet
def AirInlet_OD(SV1,SV2,AirIn,EE):

    if SV1.Ma >1:
        PD = 1-0.075*math.pow((Ma-1),1.35)
    else:
        PD = 1

    if EE.H>=0 and EE.H <=11000:
        SV1.P = 101325*math.pow((1-EE.H/44.331),5.25588)
        SV1.T = EE.T0 - 6.5*EE.H
    else:
        SV1.P = 22632*math.pow(e,(11-H)/6342)
        SV1.T = 216.62+EE.T0-288.15

    SV2.P = SV1.P*PD
    SV2.T = SV1.T

def Fan_OD(SV1,SV2,Fan):
    global frac_Comp
    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
    fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
    frac_Comp = fluid.get_mole_fractions()
    #Comp_TRAT_OD =  math.sqrt(T_in / TSL)
    Fan.CN_OD = (Fan.N  / math.sqrt(SV1.T/TSL))/Fan.N_ref
    PRCSF_local = Fan.PR_SFDP * (1.0 - DEDG.DGF)
    WACSF_local = Fan.CM_SFDP * (1.0 - DEDG.DGF)
    ETACSF_local = Fan.EF_SFDP * (1.0 - DEDG.DEF)

    if Fan.Adap:
        DeltaCN = (Fan.CN_DP - Fan.CN_OD) / Fan.CN_DP
        Fan.PR_SFOD = 1.0 + Fan.PR_b * DeltaCN + Fan.PR_c * DeltaCN * DeltaCN
        Fan.CM_SFOD = 1.0 + Fan.CM_b * DeltaCN + Fan.CM_c * DeltaCN * DeltaCN
        Fan.EF_SFOD = 1.0 + Fan.EF_b * DeltaCN + Fan.EF_c * DeltaCN * DeltaCN

        PRCSF_local *= Fan.PR_SFOD 
        WACSF_local *= Fan.CM_SFOD
        ETACSF_local *= Fan.EF_SFOD


    SV2.P = SV1.P*Fan.PR_OD
    Fan.PR_AC = (Fan.PR_OD - 1.0) / PRCSF_local + 1.0
    (Fan.PR_AC, Fan.CM_AC, Fan.EF_AC) = ReadMap(Fan.CN_OD, Fan.PR_AC,0,  Fan_GTMap)

    Fan.EF_OD = ETACSF_local * Fan.EF_AC
    Fan.CM_OD = WACSF_local * Fan.CM_AC
    SV1.MF = Fan.CM_OD * (SV1.P/PSL) / math.sqrt(SV1.T/TSL)
    SV1.H = PropsSI('H','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV1.S = PropsSI('S','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_oI = PropsSI('H','P',SV2.P,'S',SV1.S,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV2.H = SV1.H + (H_oI - SV1.H) / Fan.EF_OD
    SV2.T = PropsSI('T','H',SV2.H,'P',SV2.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SW = SV2.H - SV1.H #comp specific work
    #CW = SV1cal.MF * SW #comp total work (gcl*Wc1/EFtc)
    Fan.CW = SV1.MF * SW
    SV2.MF = SV1.MF
    Fan.CW = SV2.MF*SW


def LPC_OD(SV1,SV2,LPC):
    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
    fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
    frac_Comp = fluid.get_mole_fractions()
    #Comp_TRAT_OD =  math.sqrt(T_in / TSL)
    LPC.CN_OD = (LPC.N  / math.sqrt(SV1.T/TSL))/LPC.N_ref
    PRCSF_local = LPC.PR_SFDP * (1.0 - DEDG.DGB)
    WACSF_local = LPC.CM_SFDP * (1.0 - DEDG.DGB)
    ETACSF_local = LPC.EF_SFDP * (1.0 - DEDG.DEB)
    if LPC.Adap:
        DeltaCN = (LPC.CN_DP - LPC.CN_OD) / LPC.CN_DP
        LPC.PR_SFOD = 1.0 + LPC.PR_b * DeltaCN + LPC.PR_c * DeltaCN * DeltaCN
        LPC.CM_SFOD = 1.0 + LPC.CM_b * DeltaCN + LPC.CM_c * DeltaCN * DeltaCN
        LPC.EF_SFOD = 1.0 + LPC.EF_b * DeltaCN + LPC.EF_c * DeltaCN * DeltaCN

        PRCSF_local *= LPC.PR_SFOD 
        WACSF_local *= LPC.CM_SFOD
        ETACSF_local *= LPC.EF_SFOD


    SV2.P = SV1.P*LPC.PR_OD
    LPC.PR_AC = (LPC.PR_OD - 1.0) / PRCSF_local + 1.0
    (LPC.PR_AC, LPC.CM_AC, LPC.EF_AC) = ReadMap(LPC.CN_OD, LPC.PR_AC,0,  LPC_GTMap)

    LPC.EF_OD = ETACSF_local * LPC.EF_AC
    LPC.CM_OD = WACSF_local * LPC.CM_AC
    SV1.MF = (LPC.CM_OD*(SV1.P/PSL))/(math.sqrt(SV1.T/TSL))
    SV1.H = PropsSI('H','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV1.S = PropsSI('S','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_oI = PropsSI('H','P',SV2.P,'S',SV1.S,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV2.H = SV1.H + (H_oI - SV1.H) / LPC.EF_OD
    SV2.T = PropsSI('T','P',SV2.P,'H',SV2.H,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SW = SV2.H - SV1.H #comp specific work
    #CW = SV1cal.MF * SW #comp total work (gcl*Wc1/EFtc)
    LPC.CW = SV1.MF * SW
    SV2.MF = SV1.MF
    LPC.CW = SV2.MF*SW  

def Divide_OD(SV1,SV2,SV21,Divide):
    Divide.CMCB = Divide.WACB/Divide.valve_DP*Divide.valve
    SV2.MF = SV1.MF *Divide.CMCB
    SV2.T = SV1.T
    SV2.P = SV1.P
    SV21.T = SV1.T
    SV21.P = SV1.P
    SV21.MF = SV1.MF-SV2.MF


def Nozzle(SV1,SV2,Nozzle):#Pt,Tt,

    SV2.Pt = SV1.Pt*Nozzle.Pr
    SV2.Tt = SV1.Tt
    SV2.D = PropsSI('D','P',SV2.Pt,'T',SV2.Tt,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac))
    SV2.Wc = SV2.MF/SV2.D
    SV2.Cp = PropsSI('C','P',SV2.Pt,'T',SV2.Tt,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac))
    SV2.Cv = PropsSI('O','P',SV2.Pt,'T',SV2.Tt,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac))
    SV2.H= PropsSI('H','P',SV2.Pt,'T',SV2.Tt,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac))
    SV2.K = SV2.Cp/SV2.Cv
    Pnz = SV2.Pt/Pamb
    Pnzcr = math.pow(((1+SV2.K)/2),(SV2.K/(SV2.K-1)))

    if Pnz<=Pnzcr:
        SV2.P = Pamb
        SV2.V = Nozzle.Vc*math.sqrt(2*SV2.Cp*SV2.Tt*(1-math.pow(Pnz,(1-SV2.K)/SV2.K)))#Nozzle.Vc != rrr???
        SV2.Ma = math.pow(Pnz,(SV2.K-1)/SV2.K)
        rrr = math.sqrt(((SV2.K+1)*SV2.Ma*SV2.Ma/2)/(1+(SV2.K-1)*SV2.Ma*SV2.Ma/2))#速度系数
        qqq =math.pow((SV2.K+1)/2,1/(SV2.K-1))*rrr*math.pow((1-(SV2.K-1)/(SV2.K+1)*rrr*rrr),(1/(SV2.K-1))) #无量纲密度流
        KKK = math.sqrt(SV2.K/8.314*math.pow(2/(SV2.K+1),(SV2.K+1)/(SV2.K-1)))
        SV2.MF = (KKK*SV2.Pt*Nozzle.Aout*qqq)/math.sqrt(SV2.Tt)
    
    else:
        SV2.P = SV2.Pt/Pnzcr
        SV2.V = Nozzle.Vc*math.sqrt(2*SV2.Cp*SV2.Tt*(1-math.pow(Pnzcr,(1-SV2.K)/SV2.K)))
        SV2.Ma = 1.0
        rrr = math.sqrt(((SV2.K+1)*SV2.Ma*SV2.Ma/2)/(1+(SV2.K-1)*SV2.Ma*SV2.Ma/2))#速度系数
        qqq =math.pow((SV2.K+1)/2,1/(SV2.K-1))*rrr*math.pow((1-(SV2.K-1)/(SV2.K+1)*rrr*rrr),(1/(SV2.K-1))) #无量纲密度流
        KKK = math.sqrt(SV2.K/8.314*math.pow(2/(SV2.K+1),(SV2.K+1)/(SV2.K-1)))
        SV2.MF = (KKK*SV2.Pt*Nozzle.Aout*qqq)/math.sqrt(SV2.Tt)

#CE HPC
def HPC_OD(SV1,SV2,HPC,N=1,SV11=0,SV12=0,SV13=0,SV14=0):
    global frac_Comp
    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
    fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
    frac_Comp = fluid.get_mole_fractions()
    #Comp_TRAT_OD =  math.sqrt(T_in / TSL)
    HPC.CN_OD = (HPC.N  / math.sqrt(SV1.T/TSL))/HPC.N_ref
    PRCSF_local = HPC.PR_SFDP * (1.0 - DEDG.DGC)
    WACSF_local = HPC.CM_SFDP * (1.0 - DEDG.DGC)
    #ETACSF_local = HPC.EF_SFDP * (1.0 - DEDG.DEC)
    ETACSF_local = 1 * (1.0 - DEDG.DEC)

    if HPC.Adap:
        DeltaCN = (HPC.CN_DP - HPC.CN_OD) / HPC.CN_DP
        HPC.PR_SFOD = 1.0 + HPC.PR_b * DeltaCN + HPC.PR_c * DeltaCN * DeltaCN
        HPC.CM_SFOD = 1.0 + HPC.CM_b * DeltaCN + HPC.CM_c * DeltaCN * DeltaCN
        HPC.EF_SFOD = 1.0 + HPC.EF_b * DeltaCN + HPC.EF_c * DeltaCN * DeltaCN

        PRCSF_local *= HPC.PR_SFOD 
        WACSF_local *= HPC.CM_SFOD
        ETACSF_local *= HPC.EF_SFOD


    SV2.P = SV1.P*HPC.PR_OD
    HPC.PR_AC = (HPC.PR_OD - 1.0) / PRCSF_local + 1.0
    (HPC.PR_AC, HPC.CM_AC, HPC.EF_AC) = ReadMap(HPC.CN_OD, HPC.PR_AC,0,  HPC_GTMap)

    HPC.EF_OD = ETACSF_local * HPC.EF_AC
    HPC.CM_OD = WACSF_local * HPC.CM_AC
    #SV1.MF = HPC.CM_OD * (SV1.P/PSL) / math.sqrt(SV1.T/TSL)
    SV1.MF = HPC.CM_OD * (SV1.P) / math.sqrt(SV1.T)/1000
    SV1.H = PropsSI('H','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV1.S = PropsSI('S','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    H_oI = PropsSI('H','P',SV2.P,'S',SV1.S,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV2.H = SV1.H + (H_oI - SV1.H) / HPC.EF_OD
    SV2.T = PropsSI('T','P',SV2.P,'H',SV2.H,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))

    SW = SV2.H - SV1.H #comp specific work
    HPC.CW = SV1.MF * SW
    SV2.MF = SV1.MF
    HPC.CW = SV2.MF*SW

    #CW = SV1cal.MF * SW #comp total work (gcl*Wc1/EFtc)
    # if N ==1:
    #     HPC.CW = SV1.MF * SW
    #     SV2.MF = SV1.MF
    #     HPC.CW = SV2.MF*SW
    # elif N ==2:
    #     bypass(SV1,SV2,SV11,HPC)
    #     Bleed(SV1,SV11,HPC,1)
    #     SV2.MF = SV1.MF-SV11.MF
    #     HPC.CW = SV2.MF*SW+HPC.CWB1
    # elif N == 3:
    #     bypass(SV1,SV2,SV11,HPC,3,SV12)
    #     Bleed(SV1,SV11,HPC,1)
    #     Bleed(SV1,SV12,HPC,2)
    #     SV2.MF = SV1.MF-SV11.MF-SV12.MF
    #     HPC.CW = SV2.MF*SW+HPC.CWB1+HPC.CWB2
    # elif N ==4:
    #     bypass(SV1,SV2,SV11,HPC,4,SV12,SV13)
    #     Bleed(SV1,SV11,HPC,1)
    #     Bleed(SV1,SV12,HPC,2)
    #     Bleed(SV1,SV13,HPC,3)
    #     SV2.MF = SV1.MF-SV11.MF-SV12.MF-SV13.MF
    #     HPC.CW = SV2.MF*SW+HPC.CWB1+HPC.CWB2+HPC.CWB3
    # elif N ==5:
    #     bypass(SV1,SV2,SV11,HPC,5,SV12,SV13,SV14)
    #     Bleed(SV1,SV11,HPC,1)
    #     Bleed(SV1,SV12,HPC,2)
    #     Bleed(SV1,SV13,HPC,3)
    #     Bleed(SV1,SV14,HPC,4)
    #     SV2.MF = SV1.MF-SV11.MF-SV12.MF-SV13.MF-SV14.MF
    #     HPC.CW = SV2.MF*SW+HPC.CWB1+HPC.CWB2+HPC.CWB3+HPC.CWB4
#CE Combustor
def Burn_OD(SV1,SV2,SVF,Burn):
    global frac
    Burn.PD_OD = Burn.PD_Ref * SV1.MF * SV1.MF *SV1.T / SV1.P
    Burn.PD_OD = 0.03
    SV2.P = SV1.P * (1.0 - Burn.PD_OD)

    #SV2.H = (SV1.MF* SV1.H+ SVF.MF*Burn.LHV+ Burn.FuelEnthalpy*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
    SV2.H = (SV1.MF* SV1.H+ Burn.FuelEnthalpy*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
    #CB_P.EF = (SV2.H*SV2.MF)/(SV1.H*SV1.MF+SVF.MF*CB_P.LHV+CB_P.FuelEnthalpy*SVF.MF)
    frac = GasProp(SV1.MF, SVF.MF)

    SV2.T = PropsSI('T','P',SV2.P,'H',SV2.H,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    #HT = T_ou - T_in
#        (HT, Burn_CM_R, Burn_ETA_R) = ReadMap(0.00001 * P_in, HT, 0, LM2500Map.FixedBURNERMap1)
#        errorLocal = abs(Burn_ETA_R * ETABSF_DP - Burn_OD_Eff)
#        if errorLocal < error:
#            break
#        else:
#            Burn_OD_Eff = Burn_OD_Eff - 0.6 * errorLocal   

    SV2.MF = SV1.MF + SVF.MF
#CE HPT
def HPT_OD(SV1,SV2,HPT):
    HPT.CN_OD = (HPT.N  / math.sqrt(SV1.T/TSL)) / HPT.N_ref
    ####加入取小数位后5位
    PRCSF_local = HPT.PR_SFDP * (1.0 - DEDG.DGT+0.27)
    WACSF_local = HPT.CM_SFDP * (1.0 - DEDG.DGT+0.27)
    #ETACSF_local = HPT.EF_SFDP * (1.0 - DEDG.DET)
    ETACSF_local = 1 * (1.0 - DEDG.DET-0.369)
    if HPT.Adap:
        DeltaCN = (HPT.CN_DP - HPT.CN_OD) / HPT.CN_DP
        HPT.PR_SFOD = 1.0 + HPT.PR_b * DeltaCN + HPT.PR_c * DeltaCN * DeltaCN
        HPT.CM_SFOD = 1.0 + HPT.CM_b * DeltaCN + HPT.CM_c * DeltaCN * DeltaCN
        HPT.EF_SFOD = 1.0 + HPT.EF_b * DeltaCN + HPT.EF_c * DeltaCN * DeltaCN

        PRCSF_local *= HPT.PR_SFOD
        WACSF_local *= HPT.CM_SFOD
        ETACSF_local *= HPT.EF_SFOD

    SV2.P = SV1.P/HPT.PR_OD
    HPT.PR_AC = (HPT.PR_OD - 1.0) / PRCSF_local + 1.0
    (HPT.PR_AC, HPT.CM_AC, HPT.EF_AC) = ReadMap(HPT.CN_OD, HPT.PR_AC, 0, HPT_GTMap)

    HPT.EF_OD = ETACSF_local * HPT.EF_AC
    HPT.CM_OD = WACSF_local * HPT.CM_AC

    #SV1.MF = HPT.CM_OD * (SV1.P/PSL) / math.sqrt(SV1.T/TSL)
    SV2.MF = HPT.CM_OD * (SV1.P) / math.sqrt(SV1.T)/1000
    #SV2.MF = SV1.MF
    #T_in = 1499.9464246355676#T_in 
    #P_in = 1759656.3814859465
    #frac = [0.1229331109692631, 0.7505839086387777, 0.040145246737652754, 0.07735476043568186, 0.008982973218624533]

    SV1.S = PropsSI('S','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    SV1.H = PropsSI('H','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    H_oI = PropsSI('H','P',SV2.P,'S',SV1.S,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    SV2.H = SV1.H - (SV1.H - H_oI) * HPT.EF_OD #(SV1.H - SV2.H)/(SV1.H - H_oI) = HPT_EFF_DP
    SV2.T = PropsSI('T','P',SV2.P,'H',SV2.H,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

    HPT.CW = SV1.MF * (SV1.H - SV2.H)
    #SV2.MF = SV1.MF

def LPT_OD(SV1,SV2,LPT):
    LPT.CN_OD = (LPT.N  / math.sqrt(SV1.T/TSL)) / LPT.N_ref
    PRCSF_local = LPT.PR_SFDP * (1.0 - DEDG.DGP)
    WACSF_local = LPT.CM_SFDP * (1.0 - DEDG.DGP)
    ETACSF_local = LPT.EF_SFDP * (1.0 - DEDG.DEP)

    if LPT.Adap:
        DeltaCN = (LPT.CN_DP - LPT.CN_OD) / LPT.CN_DP
        LPT.PR_SFOD = 1.0 + LPT.PR_b * DeltaCN + LPT.PR_c * DeltaCN * DeltaCN
        LPT.CM_SFOD = 1.0 + LPT.CM_b * DeltaCN + LPT.CM_c * DeltaCN * DeltaCN
        LPT.EF_SFOD = 1.0 + LPT.EF_b * DeltaCN + LPT.EF_c * DeltaCN * DeltaCN

        PRCSF_local *= LPT.PR_SFOD
        WACSF_local *= LPT.CM_SFOD
        ETACSF_local *= LPT.EF_SFOD

    SV2.P = SV1.P/LPT.PR_OD
    LPT.PR_AC = (LPT.PR_OD - 1.0) / PRCSF_local + 1.0
    (LPT.PR_AC, LPT.CM_AC, LPT.EF_AC) = ReadMap(LPT.CN_OD, LPT.PR_AC, 0, LPT_GTMap)

    LPT.EF_OD = ETACSF_local * LPT.EF_AC
    LPT.CM_OD = WACSF_local * LPT.CM_AC
    SV1.MF = LPT.CM_OD * (SV1.P/PSL) / math.sqrt(SV1.T/TSL)
    SV2.MF = SV1.MF
    #T_in = 1499.9464246355676#T_in 
    #P_in = 1759656.3814859465
    #frac = [0.1229331109692631, 0.7505839086387777, 0.040145246737652754, 0.07735476043568186, 0.008982973218624533]

    SV1.S = PropsSI('S','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    SV1.H = PropsSI('H','T',SV1.T,'P',SV1.P,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    H_oI = PropsSI('H','P',SV2.P,'S',SV1.S,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    SV2.H = SV1.H - (SV1.H - H_oI) * LPT.EF_OD #(SV1.H - SV2.H)/(SV1.H - H_oI) = LPT_EFF_DP
    SV2.T = PropsSI('T','P',SV2.P,'H',SV2.H,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

    LPT.CW = SV1.MF * (SV1.H - SV2.H)
    SV2.MF = SV1.MF

#CE Bleed Lost
def Lost_OD(SV1,SV2,SV21,Lost):
    Lost.CMCB = Lost.WACB/Lost.valve_DP*Lost.valve
    ClassAF_Assignment(SV1,SV2)
    ClassAF_Assignment(SV1,SV21)    
    SV2.MF = SV1.MF *Lost.CMCB
    SV21.MF = SV1.MF-SV2.MF

def Mix_OD(SV1,SV2,SV11,MIX_DP,SVF):#cal T and P of SV2 #ingore EF and DP
    global frac
    SV2.MF = SV1.MF+SV11.MF
    SV2.P = SV1.P*(1-MIX_DP.PD_DP)
    SV11.H = PropsSI('H','P',SV11.P,'T',SV11.T,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
    SV1.H = PropsSI('H','P',SV1.P,'T',SV1.T,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    SV2.H = (SV1.H*SV1.MF+SV11.H*SV11.MF)*MIX_DP.EF_DP/(SV2.MF)
    frac = GasProp(SV1.MF - SVF.MF + SV11.MF, SVF.MF)#update frac
    SV2.T = PropsSI('T','H',SV2.H,'P',SV2.P,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
    #H_oI = (SV1.H * SV1.MF + SV11.H * SV11.MF)/(SV1.MF+SV11.MF)
    #MIX_DP.EF = H_oI*(SV1.MF+SV11.MF)/(SV2.H*SV2.MF)
    #MIX_DP.PD = 0
    
def MixAir_OD(SV1,SV2,SV11,MIX_DP):#cal T and P of SV2 #ingore EF and DP
    SV2.MF = SV1.MF+SV11.MF
    SV2.P = SV1.P*(1-MIX_DP.PD_DP)
    SV11.H = PropsSI('H','P',SV11.P,'T',SV11.T,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
    SV1.H = PropsSI('H','P',SV1.P,'T',SV1.T,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
    SV2.H = (SV1.H*SV1.MF+SV11.H*SV11.MF)*MIX_DP.EF_DP/(SV2.MF)
    SV2.T = PropsSI('T','H',SV2.H,'P',SV2.P,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')

def FanBypass_OD(SV2,SV21,SV11,FanBypass):
    SV21.P = SV2.P
    SV21.T = SV2.T    
    SV11.P = SV2.P
    SV11.T = SV2.T
    SV21.MF = SV2.MF*1/(FanBypass.BR_DP+1)
    SV11.MF = SV2.MF*FanBypass.BR_DP/(FanBypass.BR_DP+1)

def PD_OD(SV1,SV2,PD):
    ClassAF_Assignment(SV1,SV2)
    SV2.P = (1-PD.PD_DP)*SV1.P
