# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 18:55:01 2019

@author: dell
"""
class Nozzle():#79
    def __init__(self):
        self.Ain = 0.0
        self.Aout = 0.0
        self.Vc =1.0 #speed coefficient
        self.Pr = 1.0 #total pressure recovery coefficient
        pass
    
class FanBypass():
    def __init__(self):
        self.BR_DP = 0.0#Bypass ratio
        self.BR_OD = 0.0

class PressureDamage():
    def __init__(self):
        self.PD_DP = 0.0#Bypass ratio
        self.PD_OD = 0.0#Bypass ratio        


class EE_DPIP():
    def __init__(self):
        self.Ain = 0.0
        self.H = 0.0
        self.EF_tol = 0.0
        self.EF_Fuel = 0.0
        self.Thrust=0.0
        self.Thrust_ExCul = 0.0#external culvert
        self.Thrust_InCul = 0.0#internal culvert
        self.Speed = 0.0
        self.WL = 0.0#Whight_load
        self.T0 = 0.0#tem in the zero h
        pass


class Com_DPIP():
    def __init__(self):
        self.PR_DP= 0.0
        self.EF_DP = 0.0
        self.CM_DP = 0.0
        self.CW_DP = 0.0

        self.PR_SFDP = 0.0
        self.EF_SFDP = 0.0
        self.CM_SFDP = 0.0

        self.PR_SFOD = 0.0
        self.EF_SFOD = 0.0
        self.CM_SFOD = 0.0
        #data in original line
        self.PR_OD = 0.0
        self.EF_OD = 0.0
        self.CM_OD = 0.0
        self.CN_OD = 0.0
        #data after revised
        self.PR_AC = 0.0
        self.EF_AC = 0.0
        self.CM_AC = 0.0
        self.CN_AC = 0.0

        self.SM_DP = 0.0

        self.PR_MAP = 0.0
        self.EF_MAP = 0.0
        self.CM_MAP = 0.0

        self.DE = 0.0
        self.DG = 0.0

        self.N_ref = 0.0
        self.N = 0.0
        self.CN_DP = 0.0


        self.PRB1 = 0.0
        self.WACB1 = 0.0
        self.PRRB1 = 0.0
        self.EFB1 = 0.0
        self.CWB1 = 0.0
        self.valve1_DP = 1
        self.valve1 = 0.0
        self.CMCB1 = 0.0

        self.PRB2 = 0.0
        self.WACB2 = 0.0
        self.PRRB2 = 0.0
        self.EFB2 = 0.0
        self.CWB2 = 0.0
        self.valve2_DP = 1
        self.valve2 = 0.0
        self.CMCB2 = 0.0

        self.PRB3 = 0.0
        self.WACB3 = 0.0
        self.PRRB3 = 0.0
        self.EFB3 = 0.0
        self.CWB3 = 0.0
        self.valve3_DP = 1
        self.valve3 = 0.0
        self.CMCB3 = 0.0

        self.PRB4 = 0.0
        self.WACB4 = 0.0
        self.PRRB4 = 0.0
        self.EFB4 = 0.0
        self.CWB4 = 0.0
        self.valve4_DP = 1
        self.valve4 = 0.0
        self.CMCB4 = 0.0

        self.Adap = True
        self.PR_b = 0.0
        self.PR_c = 0.0
        self.EF_b = 0.0
        self.EF_c = 0.0
        self.CM_b = 0.0
        self.CM_c = 0.0

        self.MapNo = 1

        self.CW = 0.0
        pass

class Burn_DPIP():
    def __init__(self):
        self.EF = 0.0#effi
        self.PD_DP = 0.0#pressure damage
        self.PD_OD = 0.0
        self.PD_AC = 0.0
        self.PD_Ref = 0.0
        self.LHV = 0.0
        self.FuelEnthalpy = 0.0
        self.FuelF = 0.0
        #self.c1 = 0.0
        #self.c2 = 0.0
        #self.c3 = 0.0
        #self.c4 = 0.0
        #self.ic4 = 0.0
        #self.c5 = 0.0
        #self.ic5 = 0.0
        #self.c6 = 0.0
        #self.N2 = 0.0
        #self.co2 = 0.0
        pass

class Lost_DPIP():
    def __init__(self):
        self.WACB = 0.0 #dp CMCB
        self.CMCB = 0.0
        self.WLost = 0.0
        self.valve_DP = 0.0
        self.valve = 0.0
        pass

class Turb_DPIP():
    def __init__(self):
        self.PR_DP= 0.0
        self.EF_DP = 0.0
        self.CM_DP = 0.0
        self.CW_DP = 0.0

        self.PR_SFDP = 0.0
        self.EF_SFDP = 0.0
        self.CM_SFDP = 0.0

        self.PR_SFOD = 0.0
        self.EF_SFOD = 0.0
        self.CM_SFOD = 0.0
        
        self.PR_OD = 0.0
        self.EF_OD = 0.0
        self.CM_OD = 0.0
        self.CN_OD = 0.0
        self.PR_AC = 0.0
        self.EF_AC = 0.0
        self.CM_AC = 0.0
        self.CN_AC = 0.0
        self.SM_DP = 0.0

        self.PR_MAP = 0.0
        self.EF_MAP = 0.0
        self.CM_MAP = 0.0

        self.DE = 0.0
        self.DG = 0.0

        self.N_ref = 0.0
        self.N = 0.0
        self.CN_DP = 0.0


        self.Adap = True
        self.PR_b = 0.0
        self.PR_c = 0.0
        self.EF_b = 0.0
        self.EF_c = 0.0
        self.CM_b = 0.0
        self.CM_c = 0.0
        self.P_ou = 0.0

        self.MapNo = 1

        self.TW = 0.0
        pass
    pass

class Mix_DPIP():
    def __init__(self):
        self.EF_DP = 0.0
        self.PD_DP = 0.0
        self.EF_OD = 0.0
        self.PD_OD = 0.0
        pass

class DEDG():
    def __init__(self):
        self.DGF =0 #xone[2]
        self.DEF =0 #xone[3]
        self.DGB =0 #xone[2]
        self.DEB =0 #xone[3]
        self.DGC =0 #xone[2]
        self.DEC =0 #xone[3]
        self.DGT =0 #xone[4]
        self.DET =0 #xone[5]
        self.DGP =0 #xone[6]
        self.DEP =0 #xone[7]
        pass