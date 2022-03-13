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
import csv
import xlrd
import pandas as pd
from pandas import read_csv



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
    #质量分数
    c1_m=0.91886
    c2_m=0.05576
    c3_m=0.0074
    c4_m=0
    ic4_m=0
    c5_m=0
    ic5_m=0
    c6_m=0
    N2_m=0.00292
    co2_m=0.01505
    #摩尔分数
    totalMass_F=c1_m+c2_m+c3_m+N2_m+co2_m
    c1=c1_m/totalMass_F
    c2=c2_m/totalMass_F
    c3=c3_m/totalMass_F
    c4=0
    ic4=0
    c5=0
    ic5=0
    c6=0
    N2=N2_m/totalMass_F
    co2=co2_m/totalMass_F

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
def run(input):

    book = xlrd.open_workbook(input)
    sheet1 = book.sheets()[0]
    nrows = sheet1.nrows



    fielddata_NH = sheet1.col_values(0) #转速
    fielddata_MFF = sheet1.col_values(1)#kg/s #燃料流量
    fielddata_T25 = sheet1.col_values(2) #压气机入口温度
    fielddata_P25 = sheet1.col_values(3) #压气机入口压力
    fielddata_T3 = sheet1.col_values(4)  #压气机出口温度
    fielddata_P3 = sheet1.col_values(5)  #压气机出口压力
    fielddata_T45 = sheet1.col_values(6) #透平出口温度
    fielddata_P45 = sheet1.col_values(7) #透平出口压力
    fielddata_dp=sheet1.col_values(8)  #进气压差
    fielddata_p0=sheet1.col_values(9)  #大气压力
    fielddata_IGV=sheet1.col_values(10)  #IGV开度
    fielddata_M1=sheet1.col_values(11)  #压气机流量
    fielddata_TF=sheet1.col_values(12)  #天然气温度
    fielddata_PF=sheet1.col_values(13)  #天然气压力
    # fielddata_HPCB1_valve = 1
    # fielddata_HPCB2_valve = 1
    # fielddata_HPTRotorDivide_valve = 1
    # fielddata_HPTNGVDivide_valve = 1

    dpm10=0.0002  #压差/流量理想值
   
    result = numpy.zeros((nrows,6))

    localtime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime()) 
    dataFrame=pd.DataFrame(columns=['EFF_C','EFF_C_afterwashing','MF_C','MF_C_afterwashing','EFF_ALL','EFF_ALL_afterwashing'])
    # dataFrame.to_csv('性能计算'+localtime+'.csv',)
    #numpy.savetxt('diag'+localtime+'.csv',diag,delimiter=',')

    dataFrame2=pd.DataFrame(columns=['EFF_ALL_dp','EFF_ALL_dp0','P_dp','P_dp0'])
    for i in range(3,nrows):

        N = 3000
        T1 = fielddata_T25 [i]+273.15
        P1 = fielddata_P25 [i]*100000
        T2 = fielddata_T3 [i]+273.15
        P2 = fielddata_P3 [i]*100000+101325
        T4 = fielddata_T45 [i]+273.15
        P4 = fielddata_P45 [i]*100000
        MF = fielddata_MFF[i]
        M1=fielddata_M1[i]
        IGV=fielddata_IGV[i]
        TF=fielddata_TF[i]+273.15
        PF=fielddata_PF[i]*100000
        dp=fielddata_dp[i]  #mmh2o
        runningHour=200  #运行时间，小时

        #压气机效率计算
        global frac_Comp
        fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        #fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        fluid.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        frac_Comp = fluid.get_mole_fractions()

        H1 = PropsSI('H','T',T1,'P',P1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        S1 = PropsSI('S','T',T1,'P',P1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        H2_oI = PropsSI('H','P',P2,'S',S1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        H2 = PropsSI('H','T',T2,'P',P2,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        EFF_C=(H2_oI-H1)/(H2-H1)  #压气机效率
        W_C= M1*(H2-H1)  #压气机耗功

        #燃烧室计算
        global frac
        dp_comb = 0.03
        EFF_comb=1
        P3 = P2 * (1.0 - dp_comb)
        #质量分数
        c1_m=0.91886
        c2_m=0.05576
        c3_m=0.0074
        N2_m=0.00292
        co2_m=0.01505
        #获取组分
        fuel=CoolProp.AbstractState('REFPROP','methane&ethane&propane&nitrogen&co2')
        fuel.set_mass_fractions([c1_m, c2_m, c3_m, N2_m,co2_m])
        frac_f = fuel.get_mole_fractions()

        FuelEnthalpy=PropsSI('H','P',PF,'T',TF,'REFPROP::methane[{}]&ethane[{}]&propane[{}]&nitrogen[{}]&co2[{}]'.format(*frac_f))

        fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        #fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        fluid.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        frac_Comp = fluid.get_mole_fractions()

        FF_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::methane[{}]&ethane[{}]&propane[{}]&nitrogen[{}]&co2[{}]'.format(*frac_f))
        F_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&argon[{}]'.format(*frac_Comp))
        frac = GasProp(M1, MF)
        D_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        #H2=PropsSI('H','T',1257.6+273.15,'P',15.138*100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        #LHV=((H2-D_H0)*((SV1.MF + SVF.MF))-SV1.MF* (SV1.H-F_H0)-(Burn.FuelEnthalpy-FF_H0)*SVF.MF)/SVF.MF
        LHV=52197605.4

        #SV2.H = (SV1.MF* SV1.H+ SVF.MF*Burn.LHV+ Burn.FuelEnthalpy*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
        H3 = (M1* (H2-F_H0)+ (FuelEnthalpy-FF_H0)*MF+MF*LHV) * EFF_comb/(M1 + MF) +D_H0
        #CB_P.EF = (SV2.H*SV2.MF)/(SV1.H*SV1.MF+SVF.MF*CB_P.LHV+CB_P.FuelEnthalpy*SVF.MF)
        
        T3 = PropsSI('T','P',P3,'H',H3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        #涡轮计算

        S3 = PropsSI('S','T',T3,'P',P3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        H4_oI = PropsSI('H','P',P4,'S',S3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        H4=PropsSI('H','P',P4,'T',T4,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        EFF_T=(H3-H4)/(H3-H4_oI)
        W_T=(M1+MF)*(H3-H4)

        P=W_T-W_C

        EFF_ALL=P/(MF*LHV)
        #dataFrame=pd.DataFrame(columns=['EFF_C','EFF_C_afterwashing','MF_C','MF_C_afterwashing','EFF_ALL','EFF_ALL_afterwashing'])
        dataFrame.loc[i-2]=[EFF_C,EFF_C*1.015*runningHour/200,M1*P1/T1/T1,M1*P1/T1/T1*1.02*runningHour/200,EFF_ALL,EFF_ALL*1.02*runningHour/200]
        dataFrame.to_csv('性能计算'+localtime+'.csv')

        #改用理想进气压力计算
        dp0=dpm10*M1*M1
        
        P1=P1+(dp-dp0)*9.8 #理想进气压力
        fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        #fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        fluid.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        frac_Comp = fluid.get_mole_fractions()

        H1 = PropsSI('H','T',T1,'P',P1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        S1 = PropsSI('S','T',T1,'P',P1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        H2_oI = PropsSI('H','P',P2,'S',S1,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        H2 = PropsSI('H','T',T2,'P',P2,'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]'.format(*frac_Comp))
        EFF_C2=(H2_oI-H1)/(H2-H1)  #压气机效率
        W_C= M1*(H2-H1)  #压气机耗功
        

        #燃烧室计算
        
        dp_comb = 0.03
        EFF_comb=1
        P3 = P2 * (1.0 - dp_comb)
        #质量分数
        c1_m=0.91886
        c2_m=0.05576
        c3_m=0.0074
        N2_m=0.00292
        co2_m=0.01505
        #获取组分
        fuel=CoolProp.AbstractState('REFPROP','methane&ethane&propane&nitrogen&co2')
        fuel.set_mass_fractions([c1_m, c2_m, c3_m, N2_m,co2_m])
        frac_f = fuel.get_mole_fractions()

        FuelEnthalpy=PropsSI('H','P',PF,'T',TF,'REFPROP::methane[{}]&ethane[{}]&propane[{}]&nitrogen[{}]&co2[{}]'.format(*frac_f))

        fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        #fluid.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        fluid.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        frac_Comp = fluid.get_mole_fractions()

        FF_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::methane[{}]&ethane[{}]&propane[{}]&nitrogen[{}]&co2[{}]'.format(*frac_f))
        F_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&argon[{}]'.format(*frac_Comp))
        frac = GasProp(M1, MF)
        D_H0=PropsSI('H','T',298.15,'P',100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        #H2=PropsSI('H','T',1257.6+273.15,'P',15.138*100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        #LHV=((H2-D_H0)*((SV1.MF + SVF.MF))-SV1.MF* (SV1.H-F_H0)-(Burn.FuelEnthalpy-FF_H0)*SVF.MF)/SVF.MF
        LHV=52197605.4

        #SV2.H = (SV1.MF* SV1.H+ SVF.MF*Burn.LHV+ Burn.FuelEnthalpy*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
        H3 = (M1* (H2-F_H0)+ (FuelEnthalpy-FF_H0)*MF+MF*LHV) * EFF_comb/(M1 + MF) +D_H0
        #CB_P.EF = (SV2.H*SV2.MF)/(SV1.H*SV1.MF+SVF.MF*CB_P.LHV+CB_P.FuelEnthalpy*SVF.MF)
        
        T3 = PropsSI('T','P',P3,'H',H3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        #涡轮计算

        S3 = PropsSI('S','T',T3,'P',P3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        H4_oI = PropsSI('H','P',P4,'S',S3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        H4=PropsSI('H','P',P4,'T',T4,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))

        EFF_T=(H3-H4)/(H3-H4_oI)
        W_T=(M1+MF)*(H3-H4)

        P2=W_T-W_C

        EFF_ALL2=P2/(MF*LHV)
        

        dataFrame2.loc[i-2]=[EFF_ALL,EFF_ALL2,P,P2]
        dataFrame2.to_csv('压损影响'+localtime+'.csv')


        print(i)

    

    return [EFF_C,W_C,EFF_T,W_T,EFF_ALL,P]




#运行接口：      
result=run('input3.xlsx')