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

def Class_assignment(SV,list):
    list_len = len(list)
    if list_len ==11:
        SV.P = list[0]
        SV.T = list[1]
        SV.V = list[2]
        SV.MF  = list[3]
        SV.H = list[4]
        SV.S = list[5]
        SV.Pt = list[6]
        SV.Tt = list[7]
        SV.Rou = list[8]
        SV.As = list[9]
        SV.Ref = list[10]
    else:
        print('error:Class_assignment')

    return (SV)


def GasProp(MF_in, FuelF):   
    global frac
    # c1 = 0.88259 
    # c2 = 0.052446
    # c3 = 0.0093719
    # c4 = 0.0029607
    # ic4 = 0.0023821
    # c5 = 0.000080262
    # ic5 = 0.000929935
    # c6 = 0.0026237
    # N2 = 0.027703
    # co2 = 0.018192
    c1=0.91886
    c2=0.05576
    c3=0.0074
    c4=0
    ic4=0
    c5=0
    ic5=0
    c6=0
    N2=0.00292
    co2=0.01505
    mass_O2 = MF_in * 0.2314 - FuelF *(4 * 15.9994 / 16.0426 * c1 \
				+ 7 * 15.9994/30.0694 * c2 \
				+ 10 * 15.9994/44.0962 * c3 \
                +13 * 15.9994/58.1222 * c4 \
                + 13 * 15.9994/58.1222 * ic4 \
                +16 * 15.9994/72.15 * c5 \
                + 16 * 15.9994/72.15 * ic5 \
                +19 * 15.9994/86.17 * c6)
    mass_N2	= MF_in * 0.7552 + FuelF * N2
    mass_CO2 = MF_in * 0.0005 + FuelF *( co2 \
					+ 44.0098/16.0426 * c1 \
					+ 2 * 44.0098/30.0694 * c2 \
					+ 3* 44.0098/44.0962 * c3 \
					+4 * 44.0098/58.122 * c4 \
                    + 4 * 44.0098/58.122 * ic4 \
                    +5 * 44.0098/72.15 * c5 \
                    + 5 * 44.0098/72.15 * ic5 \
                    +6 * 44.0098/86.17 * c6)
    mass_H2O = FuelF *( 2 * 18.0152/16.0426 * c1 \
					+ 3 * 18.0152/30.0694 * c2 \
					+ 4 * 18.0152/44.0962 * c3 \
                    +5 * 18.0152/58.122 * c4 \
                    + 5 * 18.0152/58.1222 * ic4 \
                    +6 * 18.0152/72.15 * c5 \
                    + 6 * 18.0152/72.15 * ic5 \
                    +7 * 18.0152/86.17 * c6)
    mass_Ar = MF_in * 0.0129

    TotalMass = mass_O2 + mass_N2 + mass_CO2 + mass_H2O + mass_Ar
    mass_O2 = mass_O2 / TotalMass
    mass_N2 = mass_N2 / TotalMass
    mass_CO2 = mass_CO2 / TotalMass
    mass_H2O = mass_H2O / TotalMass
    mass_Ar = mass_Ar / TotalMass


    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&co2&water&argon')
    fluid.set_mass_fractions([mass_O2, mass_N2, mass_CO2, mass_H2O, mass_Ar])
    frac = fluid.get_mole_fractions()
    return(frac)
    
def Loc(a_M,Val, NumT):
    Coeff = 0
    locV = 0
    for i in range(NumT):
        if a_M[i] > Val:
            locV = i;
            Coeff = (Val - a_M[i - 1]) / (a_M[i] - a_M[i - 1]);
            break;
    return Coeff, locV

def ReadMap(CN_T, PR_T,SZ, map,CM_T = 0):
    if CN_T > map.AX[0] and CN_T <map.AX[-1]:
        C_flag=0 #实际转速在map里
    else:
        C_flag=1 #实际转速超出map
        print(C_flag)
        if CN_T <= map.AX[0]:
            CN_T = map.AX[0]
        else:
            CN_T = map.AX[-1]

    (CoeffCN, locCN) = Loc(map.AX,CN_T, map.nCN)

    PR_ODCN_MAP = numpy.zeros(map.nPOINTS)
    CMF_ODCN_MAP = numpy.zeros(map.nPOINTS)
    ETA_ODCN_MAP = numpy.zeros(map.nPOINTS)

    if CoeffCN == 0:
        for j in range(map.nPOINTS):
            PR_ODCN_MAP[j] = map.BX[locCN - 1, j];
            CMF_ODCN_MAP[j] = map.CX[locCN - 1, j];
            ETA_ODCN_MAP[j] = map.DX[locCN - 1, j];
    else:
        for j in range(map.nPOINTS):
            PR_ODCN_MAP[j] = map.BX[locCN - 1, j] + CoeffCN * (map.BX[locCN, j] - map.BX[locCN - 1, j])
            CMF_ODCN_MAP[j] = map.CX[locCN - 1, j] + CoeffCN * (map.CX[locCN, j] - map.CX[locCN - 1, j])
            ETA_ODCN_MAP[j] = map.DX[locCN - 1, j] + CoeffCN * (map.DX[locCN, j] - map.DX[locCN - 1, j])

    if SZ !=0 and PR_T == 0 and CM_T ==0:
        PR_T = PR_ODCN_MAP[0] + SZ * (PR_ODCN_MAP[-1] - PR_ODCN_MAP[0])
    if CM_T !=0:
        (CoeffCM, locCM) = Loc(CMF_ODCN_MAP,CM_T, map.nPOINTS)
        if CoeffCM == 0:
            PR_ODCN = PR_ODCN_MAP[locCM - 1];
            CMF_ODCN = CMF_ODCN_MAP[locCM - 1];
            ETA_ODCN = ETA_ODCN_MAP[locCM - 1];
        else:
            PR_ODCN = PR_ODCN_MAP[locCM - 1] + CoeffCM * (PR_ODCN_MAP[locCM] - PR_ODCN_MAP[locCM - 1])
            CMF_ODCN = CMF_ODCN_MAP[locCM - 1] + CoeffCM * (CMF_ODCN_MAP[locCM] - CMF_ODCN_MAP[locCM - 1])
            ETA_ODCN = ETA_ODCN_MAP[locCM - 1] + CoeffCM * (ETA_ODCN_MAP[locCM] - ETA_ODCN_MAP[locCM - 1])

    else:
        (CoeffPR, locPR) = Loc(PR_ODCN_MAP,PR_T, map.nPOINTS)
        if CoeffPR == 0:
            PR_ODCN = PR_ODCN_MAP[locPR - 1];
            CMF_ODCN = CMF_ODCN_MAP[locPR - 1];
            ETA_ODCN = ETA_ODCN_MAP[locPR - 1];
        else:
            PR_ODCN = PR_ODCN_MAP[locPR - 1] + CoeffPR * (PR_ODCN_MAP[locPR] - PR_ODCN_MAP[locPR - 1])
            CMF_ODCN = CMF_ODCN_MAP[locPR - 1] + CoeffPR * (CMF_ODCN_MAP[locPR] - CMF_ODCN_MAP[locPR - 1])
            ETA_ODCN = ETA_ODCN_MAP[locPR - 1] + CoeffPR * (ETA_ODCN_MAP[locPR] - ETA_ODCN_MAP[locPR - 1])
    return PR_ODCN, CMF_ODCN, ETA_ODCN

def ScaleMap(CN_T, map):
    if CN_T > map.AX[0] and CN_T <map.AX[-1]:
        C_flag=0 #实际转速在map里
    else:
        C_flag=1 #实际转速超出map
        print(C_flag)
        if CN_T <= map.AX[0]:
            CN_T = map.AX[0]
        else:
            CN_T = map.AX[-1]

    (CoeffCN, locCN) = Loc(map.AX,CN_T, map.nCN)

    PR_ODCN_MAP = numpy.zeros(map.nPOINTS)
    CMF_ODCN_MAP = numpy.zeros(map.nPOINTS)
    ETA_ODCN_MAP = numpy.zeros(map.nPOINTS)

    if CoeffCN == 0:
        for j in range(map.nPOINTS):
            PR_ODCN_MAP[j] = map.BX[locCN - 1, j];
            CMF_ODCN_MAP[j] = map.CX[locCN - 1, j];
            ETA_ODCN_MAP[j] = map.DX[locCN - 1, j];
    else:
        for j in range(map.nPOINTS):
            PR_ODCN_MAP[j] = map.BX[locCN - 1, j] + CoeffCN * (map.BX[locCN, j] - map.BX[locCN - 1, j])
            CMF_ODCN_MAP[j] = map.CX[locCN - 1, j] + CoeffCN * (map.CX[locCN, j] - map.CX[locCN - 1, j])
            ETA_ODCN_MAP[j] = map.DX[locCN - 1, j] + CoeffCN * (map.DX[locCN, j] - map.DX[locCN - 1, j])
    m = 0
    ETA_ODCN = 0
    for k in range(map.nPOINTS-1):

        if ETA_ODCN_MAP[k+1]>ETA_ODCN_MAP[k]:
            ETA_ODCN = ETA_ODCN_MAP[k+1]
            m = k+1
        else:
            pass
    PR_ODCN = PR_ODCN_MAP[m]
    CMF_ODCN = CMF_ODCN_MAP[m]


    return PR_ODCN, CMF_ODCN, ETA_ODCN

def valve(Valve,N=2):#N<=5
    if N==2:
        Valve.CMCB1 = Valve.WACB1/Valve.valve1_DP*Valve.valve1
    elif N ==3:
        Valve.CMCB1 = Valve.WACB1/Valve.valve1_DP*Valve.valve1
        Valve.CMCB2 = Valve.WACB2/Valve.valve2_DP*Valve.valve2
    elif N ==4:
        Valve.CMCB1 = Valve.WACB1/Valve.valve1_DP*Valve.valve1
        Valve.CMCB2 = Valve.WACB2/Valve.valve2_DP*Valve.valve2
        Valve.CMCB3 = Valve.WACB3/Valve.valve3_DP*Valve.valve3
    else:
        Valve.CMCB1 = Valve.WACB1/Valve.valve1_DP*Valve.valve1
        Valve.CMCB2 = Valve.WACB2/Valve.valve2_DP*Valve.valve2
        Valve.CMCB3 = Valve.WACB3/Valve.valve3_DP*Valve.valve3
        Valve.CMCB4 = Valve.WACB4/Valve.valve4_DP*Valve.valve4

def bypass(SV1,SV2,SV11,Bypass,N=2,SV12=0,SV13=0,SV14=0):#N<=5
    if N==2:
        valve(Bypass)
        SV11.MF = SV1.MF*Bypass.CMCB1
        SV2.MF = SV1.MF-SV11.MF
    elif N ==3:
        valve(Bypass,N=3)
        SV11.MF = SV1.MF*Bypass.CMCB1
        SV12.MF = SV1.MF*Bypass.CMCB2
        SV2.MF = SV1.MF-SV11.MF-SV12.MF
    elif N ==4:
        valve(Bypass,N=4)
        SV11.MF = SV1.MF*Bypass.CMCB1
        SV12.MF = SV1.MF*Bypass.CMCB2
        SV13.MF = SV1.MF*Bypass.CMCB3
        SV2.MF = SV1.MF-SV11.MF-SV12.MF-SV13.MF
    else:
        valve(Bypass,N=5)
        SV11.MF = SV1.MF*Bypass.CMCB1
        SV12.MF = SV1.MF*Bypass.CMCB2
        SV13.MF = SV1.MF*Bypass.CMCB3
        SV14.MF = SV1.MF*Bypass.CMCB4
        SV2.MF = SV1.MF-SV11.MF-SV12.MF-SV13.MF-SV14.MF
    return()

def Bleed(SV1,SV11,HPC,Bleed):
    if Bleed ==1:
        PR = HPC.PRRB1 * HPC.PR_OD
        SV11.P= SV1.P * PR
        #MF_B1 = Comp_WACB1 * SV1.MF
        HPC.EFB1 = math.pow(HPC.EF_OD,HPC.PRRB1)

        H_1I = PropsSI('H','P',SV11.P,'S',SV1.S,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
        SV11.H = SV1.H + (H_1I - SV1.H) / HPC.EFB1
        SV11.T = PropsSI('T','P',SV11.P,'H',SV11.H,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')   
        HPC.CWB1 = (SV11.H-SV1.H) * SV11.MF 
    elif Bleed ==2:
        PR = HPC.PRRB2 * HPC.PR_OD
        SV11.P= SV1.P * PR
        #MF_B2 = Comp_WACB2 * SV1.MF
        HPC.EFB2 = math.pow(HPC.EF_OD,HPC.PRRB2)

        H_1I = PropsSI('H','P',SV11.P,'S',SV1.S,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
        SV11.H = SV1.H + (H_1I - SV1.H) / HPC.EFB2
        SV11.T = PropsSI('T','P',SV11.P,'H',SV11.H,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')   
        HPC.CWB2 = (SV11.H-SV1.H) * SV11.MF 
    elif Bleed ==3:
        PR = HPC.PRRB3 * HPC.PR_OD
        SV11.P= SV1.P * PR
        #MF_B3 = Comp_WACB3 * SV1.MF
        HPC.EFB3 = math.pow(HPC.EF_OD,HPC.PRRB3)

        H_1I = PropsSI('H','P',SV11.P,'S',SV1.S,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
        SV11.H = SV1.H + (H_1I - SV1.H) / HPC.EFB3
        SV11.T = PropsSI('T','P',SV11.P,'H',SV11.H,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')   
        HPC.CWB3 = (SV11.H-SV1.H) * SV11.MF 
    elif Bleed ==4:
        PR = HPC.PRRB4 * HPC.PR_OD
        SV11.P= SV1.P * PR
        #MF_B4 = Comp_WACB4 * SV1.MF
        HPC.EFB4 = math.pow(HPC.EF_OD,HPC.PRRB4)

        H_1I = PropsSI('H','P',SV11.P,'S',SV1.S,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')
        SV11.H = SV1.H + (H_1I - SV1.H) / HPC.EFB4
        SV11.T = PropsSI('T','P',SV11.P,'H',SV11.H,'REFPROP::oxygen[0.2314]&nitrogen[0.7552]&argon[0.0129]&co2[0.0005]')   
        HPC.CWB4 = (SV11.H-SV1.H) * SV11.MF 
    else:
        pass

def ClassAF_Assignment(SV1,SV2):
        SV2.P = SV1.P#	 pressure
        SV2.T = SV1.T#	 temperature
        SV2.V = SV1.V# speed
        SV2.MF  = SV1.MF #mass flow	
        SV2.Wc = SV1.Wc #volumn flow
        SV2.H = SV1.H #enthalpy
        SV2.S = SV1.S #entropy
        SV2.Pt = SV1.Pt #total Pressure
        SV2.Tt = SV1.Tt#total temperature
        SV2.D = SV1.D #density
        SV2.As = SV1.As #sectional area
        SV2.Ma = SV1.Ma#MAa
        SV2.Cp = SV1.Cp
        SV2.Cv = SV1.Cv
        SV2.Dmass = SV1.Dmass
        SV2.Dmolar = SV1.Dmolar
        SV2.Rg = SV1.Rg
        SV2.K = SV1.K #heat capacity ratio
        SV2.Far = SV1.Far #fuel/air ratio
        SV2.War = SV1.War # water/air ratio
        SV2.Ref = SV1.Ref#Ref    
    
