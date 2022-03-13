'''
!/usr/bin/env python3.8 #Visual Studio C++ 2019不支持python3.9，安装CoolProp库报错
-*- coding: utf-8 -*-
Author: Ian
Created on Thu Nov 21 18:55:01 2019
Revised on Tue Feb 22 19:22:01 2022
'''

import numpy as np# 科学计算库，支持数组和矩阵运算，比python的嵌套列表结构更高效
import math # 数学函数库
import CoolProp # 物性参数库
# import scipy # 基于Numpy的科学计算库，可以处理插值、积分、优化、图像处理、常微分方程数值解的求解、信号处理等问题。
import time # 时间库，表述的日期范围被限定在1970-2038之间
# import datetime # 高级时间库，基于time进行了封装，提供了更多实用的函数，datetime、timedata、time、date、tzinfo
# import csv # 读写csv文件，pandas.read_csv()，pandas.to_csv()
import xlrd # 读写Excel文件，xlrd==1.2.0，2.0版本打开xlxs文件报错，pandas.read_excel()，pandas.to_excel()
import pandas as pd # 基于Numpy的复杂数据处理分析库，常用于人工智能和机器学习
import CoolProp.CoolProp as CP # 子库
# from pandas import read_csv
# from scipy import interpolate
# from scipy.optimize import fsolve
# from scipy.optimize import newton
# from scipy.optimize import fmin
from CoolProp.CoolProp import PropsSI # PropsSI物性函数
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH,'C:\\Program Files (x86)\\REFPROP\\')
# CP.get_global_param_string("REFPROP_version") # 获取REFPROP版本号

# 燃料天然气各组分摩尔质量
M_c1 = PropsSI("M","METHANE") # R50摩尔质量 kg/mol
M_c2 = PropsSI("M","ETHANE") # R170摩尔质量 kg/mol
M_c3 = PropsSI("M","PROPANE") # R290摩尔质量 kg/mol
M_c4 = PropsSI("M","BUTANE") # R600摩尔质量 kg/mol
M_ic4 = PropsSI("M","ISOBUTANE") # R600A摩尔质量 kg/mol
M_c5 = PropsSI("M","PENTANE") # R601摩尔质量 kg/mol
M_ic5 = PropsSI("M","ISOPENTANE") # R601A摩尔质量 kg/mol
M_c6 = PropsSI("M","HEXANE") # 己烷摩尔质量 kg/mol
M_ic6 = PropsSI("M","ISOHEXANE") # 异己烷摩尔质量 kg/mol
M_n2 = PropsSI("M","NITROGEN") # 氮气摩尔质量 kg/mol
M_co2 = PropsSI("M","CARBONDIOXIDE") # 二氧化碳摩尔质量 kg/mol
M_o2 = PropsSI("M","OXYGEN") # 氧气摩尔质量 kg/mol
M_h2o = PropsSI("M","WATER") # 水摩尔质量 kg/mol
M_ar = PropsSI("M","Ar") # 氩气摩尔质量 kg/mol

# 燃料天然气各组分低位热值 单位：J/kg
LHV_c1 = 50016000 
LHV_c2 = 47524800
LHV_c3 = 46340900
LHV_c4 = 45726800
LHV_ic4 = 45559400
LHV_c5 = 45352300
LHV_ic5 = 45254700
LHV_c6s = 45015100
LHV_n2 = 0
LHV_o2 = 0
LHV_co2 = 0

# 自定义烟气组分函数
def GasProp(p_atm, t_db, RH, MF_air, MF_fuel, c1, c2, c3, c4, ic4, c5, ic5, c6s, n2, o2, co2): # MF_air压气机进气质量流量，MF_fuel燃料天然气质量流量 
    global molefrac_gas # 烟气组分摩尔分数，包括氮气、氧气、二氧化碳、氩气、水蒸气
    global massfrac_air 
    '''
    干空气组分摩尔分数取自NASA报告1311（Gordon 1982）
    氮气：78.0840%
    氧气：20.9476%
    氩气：0.9365%
    二氧化碳：0.0319%
    合计：100.000%
    燃机进口干空气组分经大气湿度和大气压力修正后，得到进口空气各组分的质量分数和摩尔流量。湿度依据美国供暖、制冷和空调工程师协会的基础数据手册来计算。
    '''
    if -100 <= t_db < 0:
        k1 = -1.0214165 * pow(10, 4)
        k2 = -9.8702328
        k3 = -5.3765794 * pow(10, -3)
        k4 = -1.9202377 * pow(10, -7)
        k5 = -3.5575832 * pow(10, -10)
        k6 = -9.0344688 * pow(10, -14)
        k7 = -4.1635019
    elif 0 <= t_db <= 200:
        k1 = -1.0440397 * pow(10, 4)
        k2 = -16.27164
        k3 = -2.7022355 * pow(10, -2)
        k4 = 1.289036   * pow(10, -5)
        k5 = -2.4780681 * pow(10, -9)
        k6 = 0
        k7 = 6.5459673

    t_db = 1.8 * t_db + 491.67 
    p_v = pow(10, 6) * math.exp(k1 / t_db + k2 +k3 * t_db + k4 * pow(t_db, 2) + k5 * pow(t_db, 3) + k6 * pow(t_db, 4) + k7 * math.log(t_db, math.e))
    p_w = p_v * RH / 100
    x_DA = (p_atm - p_w) / p_atm

    # 湿空气各组分摩尔分数
    x_n2 = 0.780840 * x_DA
    x_o2 = 0.209476 * x_DA
    x_ar = 0.009365 * x_DA
    x_co2 = 0.000319 * x_DA
    x_h2o = 1 - x_DA

    # 湿空气各组分质量分数
    M_air = M_n2 * x_n2 + M_o2 * x_o2 + M_co2 * x_co2 + M_h2o * x_h2o + M_ar * x_ar # 湿空气摩尔质量
    m_n2 = M_n2 * x_n2 / M_air
    m_o2 = M_o2 * x_o2 / M_air
    m_co2 = M_co2 * x_co2 / M_air
    m_h2o = M_h2o * x_h2o / M_air
    m_ar = M_ar * x_ar / M_air
    massfrac_air = [m_o2, m_n2, m_ar, m_co2, m_h2o]

    # 烟气各组分质量流量
    M_fuel = M_c1 * c1 + M_c2 * c2 + M_c3 * c3 + M_c4 * c4 + M_ic4 * ic4 + M_c5 * c5 + M_ic5 * ic5 + M_c6 * c6s + M_n2 * n2 + M_o2 * o2 + M_co2 * co2
    mass_o2 = MF_air * m_o2 + MF_fuel * M_o2 / M_fuel * o2 - MF_fuel * (4 * M_o2 / 2 / M_c1 * c1 \
				+  7 * M_o2 / 2 / M_c2 * c2 \
				+ 10 * M_o2 / 2 / M_c3 * c3 \
                + 13 * M_o2 / 2 / M_c4 * c4 \
                + 13 * M_o2 / 2 / M_ic4 * ic4 \
                + 16 * M_o2 / 2 / M_c5 * c5 \
                + 16 * M_o2 / 2 / M_ic5 * ic5 \
                + 19 * M_o2 / 2 / M_c6 * c6s)
    mass_n2	= MF_air * m_n2 + MF_fuel * M_n2 / M_fuel * n2
    mass_co2 = MF_air * m_co2 + MF_fuel * (M_co2 / M_fuel * co2 \
					+     M_co2 / M_c1 * c1 \
					+ 2 * M_co2 / M_c2 * c2 \
					+ 3 * M_co2 / M_c3 * c3 \
					+ 4 * M_co2 / M_c4 * c4 \
                    + 4 * M_co2 / M_ic4 * ic4 \
                    + 5 * M_co2 / M_c5 * c5 \
                    + 5 * M_co2 / M_ic5 * ic5 \
                    + 6 * M_co2 / M_c6 * c6s)
    mass_h2o = MF_air * m_h2o + MF_fuel * (2 * M_h2o / M_c1 * c1 \
					+ 3 * M_h2o / M_c2 * c2 \
					+ 4 * M_h2o / M_c3 * c3 \
                    + 5 * M_h2o / M_c4 * c4 \
                    + 5 * M_h2o / M_ic4 * ic4 \
                    + 6 * M_h2o / M_c5 * c5 \
                    + 6 * M_h2o / M_ic5 * ic5 \
                    + 7 * M_h2o / M_c6 * c6s)
    mass_ar = MF_air * m_ar

    # 烟气各组分质量分数
    TotalMass = mass_o2 + mass_n2 + mass_co2 + mass_h2o + mass_ar
    mass_o2 = mass_o2 / TotalMass
    mass_n2 = mass_n2 / TotalMass
    mass_co2 = mass_co2 / TotalMass
    mass_h2o = mass_h2o / TotalMass
    mass_ar = mass_ar / TotalMass

    gas = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&co2&water&argon') # 自定义混合物烟气
    gas.set_mass_fractions([mass_o2, mass_n2, mass_co2, mass_h2o, mass_ar]) # 给定烟气各组分质量分数
    molefrac_gas = gas.get_mole_fractions() # 获取烟气各组分摩尔分数
    return [massfrac_air, molefrac_gas] # 返回值：空气质量分数，烟气摩尔分数

def run(input):
    book = xlrd.open_workbook(input) # 读取工作簿
    sheet1 = book.sheets()[0] # 选定工作表 sheet1 = book.sheet_by_index(0) 或 sheet1 = book.sheet_by_name('Sheet2')
    nrows = sheet1.nrows # 获取表的行数和列数 ncols = sheet1.ncols

    FieldData_NH = sheet1.col_values(0) # 转速, 获取第1列的值
    FieldData_MFF = sheet1.col_values(1) # 燃料流量 kg/s 
    FieldData_T25 = sheet1.col_values(2) # 压气机入口温度
    FieldData_P25 = sheet1.col_values(3) # 压气机入口压力
    FieldData_T3 = sheet1.col_values(4)  # 压气机出口温度
    FieldData_P3 = sheet1.col_values(5)  # 压气机出口压力
    FieldData_T45 = sheet1.col_values(6) # 透平出口温度
    FieldData_P45 = sheet1.col_values(7) # 透平出口压力
    FieldData_dp=sheet1.col_values(8)  # 进气压差
    FieldData_p0=sheet1.col_values(9)  # 大气压力
    FieldData_IGV=sheet1.col_values(10)  # IGV开度
    FieldData_M1=sheet1.col_values(11)  # 压气机流量
    FieldData_TF=sheet1.col_values(12)  # 天然气温度
    FieldData_PF=sheet1.col_values(13)  # 天然气压力
    FieldData_runningTime=sheet1.col_values(14) # 运行时间
    FieldData_c1 = sheet1.col_values(15) # 甲烷 摩尔分数
    FieldData_c2 = sheet1.col_values(16) # 乙烷 摩尔分数   
    FieldData_c3 = sheet1.col_values(17) # 丙烷 摩尔分数
    FieldData_c4 = sheet1.col_values(18) # 正丁烷
    FieldData_ic4 = sheet1.col_values(19) # 异丁烷
    FieldData_c5 = sheet1.col_values(20) # 正戊烷
    FieldData_ic5 = sheet1.col_values(21) # 异戊烷
    FieldData_c6s = sheet1.col_values(22) # 己烷
    FieldData_n2 = sheet1.col_values(23) # 氮气
    FieldData_o2 = sheet1.col_values(24) # 氧气
    FieldData_co2 = sheet1.col_values(25) # 二氧化碳
    FieldData_patm = sheet1.col_values(26) # 大气压力
    FieldData_tdb = sheet1.col_values(27) # 干球温度
    FieldData_RH = sheet1.col_values(28) # 相对湿度

    # FieldData_HPCB1_valve = 1 高压压气机动叶
    # FieldData_HPCB2_valve = 1
    # FieldData_HPTRotorDivide_valve = 1 高压透平转子
    # FieldData_HPTNGVDivide_valve = 1 高压透平喷嘴导叶

    k_dp0=0.002  # 压差-流量系数 ？？？
   
    result = np.zeros((nrows, 6)) # 结果保存在多维数组中

    localtime = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime()) # 获取本地时间
    dataFrame = pd.DataFrame(columns = ['EFF_comp', 'EFF_comp_afterwashing', 'M1_cor', 'M1_cor_afterwashing', 'EFF_GT', 'EFF_GT_afterwashing']) # 水洗前后性能参数：创建一个6列表格型数据结构，给出列标签
    # dataFrame.to_csv('性能计算'+localtime+'.csv',)
    # numpy.savetxt('diag'+localtime+'.csv', diag, delimiter=',') # 将array保存到txt文件；numpy.loadtxt() 将数据读出为array类型
    dataFrame0 = pd.DataFrame(columns = ['EFF_GT_dp','EFF_GT_dp0','W_GT_dp','W_GT_dp0']) # 更换进气滤芯前后燃机效率和出力对比

    # for i in range(3,nrows):
    for i in range(3,4):
        # 单位换算，转换为国际单位
        N = 3000
        T1 = FieldData_T25[i]+273.15 # 压气机进气温度 ℃ → K
        P1 = FieldData_P25[i]*100000 # 压气机进气压力 bar → Pa
        T2 = FieldData_T3[i]+273.15  
        P2 = FieldData_P3[i]*100000+101325
        T4 = FieldData_T45[i]+273.15
        P4 = FieldData_P45[i]*100000
        MF = FieldData_MFF[i]
        M1 = FieldData_M1[i]
        IGV = FieldData_IGV[i]
        TF = FieldData_TF[i] + 273.15
        PF = FieldData_PF[i] * 100000
        dp_filter = FieldData_dp[i] * 9.80665  # mmh2o → Pa 
        runningHour = FieldData_runningTime[i]  # 运行小时数 h
        c1 = FieldData_c1[i]
        c2 = FieldData_c2[i]    
        c3 = FieldData_c3[i]
        c4 = FieldData_c4[i]
        ic4 = FieldData_ic4[i]
        c5 = FieldData_c5[i]
        ic5 = FieldData_ic5[i]
        c6s = FieldData_c6s[i]
        n2 = FieldData_n2[i]
        o2 = FieldData_o2[i]
        co2 = FieldData_co2[i]
        patm = FieldData_patm[i]
        tdb = FieldData_tdb[i]
        RH = FieldData_RH[i]

        list = GasProp(patm, tdb, RH, M1, MF, c1, c2, c3, c4, ic4, c5, ic5, c6s, n2, o2, co2)
        massfrac_air = list[0]
        molefrac_gas = list[1]

        #压气机效率计算
        global molefrac_air
        air = CoolProp.AbstractState('REFPROP', 'oxygen&nitrogen&argon&co2&water')
        # air.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        # air.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        air.set_mass_fractions(massfrac_air)
        molefrac_air = air.get_mole_fractions()

        P1 = patm - dp_filter
        H1 = PropsSI('H', 'T', T1, 'P', P1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air)) # 压气机进气焓
        S1 = PropsSI('S', 'T', T1, 'P', P1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air)) # 压气机进气熵
        H2_s = PropsSI('H', 'P', P2, 'S', S1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air)) # 压气机等熵压缩排气焓
        H2 = PropsSI('H', 'T', T2, 'P', P2, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air)) # 压气机排气焓
        EFF_comp = (H2_s - H1) / (H2 - H1)  # 压气机等熵效率
        W_comp = M1 * (H2 - H1)  # 压气机耗功

        #燃烧室计算
        # global molefrac_gas
        dp_comb = 0.03 # 燃烧室压损系数 ？？？
        EFF_comb = 0.99 # 燃烧室效率 ？？？
        P3 = P2 * (1.0 - dp_comb) # 燃烧室出口压力 == 透平入口压力 
        #质量分数
        # c1_m=0.91886
        # c2_m=0.05576
        # c3_m=0.0074
        # N2_m=0.00292
        # co2_m=0.01505
        #获取组分
        # fuel=CoolProp.AbstractState('REFPROP','methane&ethane&propane&BUTANE&ISOBUTANE&PENTANE&ISOPENTANE&HEXANE&nitrogen&OXYGEN&co2')
        # fuel.set_mass_fractions([c1_m, c2_m, c3_m, N2_m,co2_m])
        # molefrac_fuel = fuel.get_mole_fractions()
        molefrac_fuel = [c1, c2, c3, c4, ic4, c5, ic5, c6s, n2, co2] # 燃料天然气各组分摩尔分数
        H_fuel = PropsSI('H', 'P', PF, 'T', TF, 'REFPROP::methane[{}]&ethane[{}]&propane[{}]&BUTANE[{}]&ISOBUTAN[{}]&PENTANE[{}]&IPENTANE[{}]&HEXANE[{}]&nitrogen[{}]&co2[{}]'.format(*molefrac_fuel)) # 燃料天然气焓

        # air = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        # #air.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        # air.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        # molefrac_air = air.get_mole_fractions()

        H0_fuel = PropsSI('H', 'T', 298.15, 'P', 101325, 'REFPROP::methane[{}]&ethane[{}]&propane[{}]&BUTANE[{}]&ISOBUTAN[{}]&PENTANE[{}]&IPENTANE[{}]&HEXANE[{}]&nitrogen[{}]&co2[{}]'.format(*molefrac_fuel)) # 燃料基准焓：15℃，101325Pa
        H0_air = PropsSI('H', 'T', 298.15, 'P', 101325, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air)) # 空气基准焓
        # frac = GasProp(patm, tdb, RH, M1, MF, c1, c2, c3, c4, ic4, c5, ic5, c6s, n2, o2, co2)
        H0_gas = PropsSI('H', 'T', 298.15, 'P', 101325, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas)) # 烟气基准焓

        #H2=PropsSI('H','T',1257.6+273.15,'P',15.138*100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        #LHV_fuel=((H2-H0_gas)*((SV1.MF + SVF.MF))-SV1.MF* (SV1.H-H0_air)-(Burn.H_fuel-H0_fuel)*SVF.MF)/SVF.MF
        M_fuel = c1 * M_c1 + c2 * M_c2 + c3 * M_c3 + c4 * M_c4 + ic4 * M_ic4 + c5 * M_c5 + ic5 * M_ic5 + c6s * M_c6 + n2 * M_n2 + o2 * M_o2 + co2 * M_co2 # 燃料天然气摩尔质量
        LHV_fuel = (c1 * M_c1 * LHV_c1 + c2 * M_c2 * LHV_c2 + c3 * M_c3 * LHV_c3 + c4 * M_c4 * LHV_c4 + ic4 * M_ic4 * LHV_ic4 + c5 * M_c5 * LHV_c5 + ic5 * M_ic5 * LHV_ic5 + c6s * M_c6 * LHV_c6s + n2 * M_n2 * LHV_n2 + o2 * M_o2 * LHV_o2 + co2 * M_co2 * LHV_co2) / M_fuel # 燃料天然气低位热值
        # LHV_fuel=52197605.4
        
        # 透平计算
        #SV2.H = (SV1.MF* SV1.H+ SVF.MF*Burn.LHV_fuel+ Burn.H_fuel*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
        H3 = (M1 * (H2 - H0_air) + MF * (H_fuel - H0_fuel) + MF * LHV_fuel) * EFF_comb / (M1 + MF) + H0_gas # 燃烧室出口烟气焓 == 透平入口烟气焓
        #CB_P.EF = (SV2.H*SV2.MF)/(SV1.H*SV1.MF+SVF.MF*CB_P.LHV_fuel+CB_P.H_fuel*SVF.MF)
        T3 = PropsSI('T', 'P', P3, 'H', H3, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas)) # 燃烧室出口烟气温度 == 透平入口烟气温度
        S3 = PropsSI('S', 'T', T3, 'P', P3, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas)) # 燃烧室出口烟气熵 == 透平入口烟气熵
        H4_is = PropsSI('H', 'P', P4, 'S', S3, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas)) # 等熵膨胀透平出口烟气焓
        H4 = PropsSI('H', 'P', P4, 'T', T4, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas)) # 透平出口烟气焓
        EFF_turb = (H3 - H4) / (H3 - H4_is) # 透平等熵效率
        # EFF_turb = (1870000 - H4) / (1870000 - H4_is) # 透平等熵效率

        W_turb = (M1 + MF) * (H3 - H4) # 透平功率

        # 燃机整体性能计算
        W_GT = W_turb - W_comp # 燃机功率
        EFF_GT = W_GT / (MF * LHV_fuel) # 燃机效率
        #dataFrame=pd.DataFrame(columns=['EFF_C','EFF_C_afterwashing','MF_C','MF_C_afterwashing','EFF_GT','EFF_ALL_afterwashing'])
        # 压气机水洗前后性能计算：将水洗后性能指标简化为运行小时数的线性函数，准确性值得商榷？？？ 折合流量公式错误？
        dataFrame.loc[i-2] = [EFF_comp, EFF_comp + 0.02 * runningHour / 500, M1 * math.sqrt(T1 / 288.15) *101325 / P1, M1 * math.sqrt(T1 / 288.15) *101325 / P1 + 0.04 * runningHour / 500, EFF_GT, EFF_GT + 0.02 * runningHour / 500] # 将水洗前后压气机效率、空气折合流量、燃机效率变化6个值输出到dataframe表格型数据结构。df.loc[a, b] loc行列索引。折合流量 m_cor=m*sqrt(t1)/p1 , 折合转速 n_cor=
        dataFrame.to_csv('性能计算2 ' + localtime + '.csv') # 将dataframe中数据输出为csv格式文件
        print(H1,H2,H3,H4,H4_is,EFF_turb,W_turb,H0_air,H0_fuel,H0_gas,H_fuel,LHV_fuel)
        # STEP2: 理想进气压损下燃机性能计算 
               
        dp0_filter = k_dp0 * M1**2 # 理想进气压损与空气流量平方成正比
        P1 = patm - dp0_filter # 压气机理想进气压力
        # air = CoolProp.AbstractState('REFPROP', 'oxygen&nitrogen&argon&co2&water')
        # #air.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        # air.set_mass_fractions(massfrac_air)
        # molefrac_air = air.get_mole_fractions()

        H1 = PropsSI('H', 'T', T1, 'P', P1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air))
        S1 = PropsSI('S', 'T', T1, 'P', P1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air))
        H2_s = PropsSI('H', 'P', P2, 'S', S1, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air))
        # H2 = PropsSI('H', 'T', T2, 'P', P2, 'REFPROP::oxygen[{}]&nitrogen[{}]&argon[{}]&co2[{}]&water[{}]'.format(*molefrac_air))
        EFF_comp0 = (H2_s - H1) / (H2 - H1)  # 压气机等熵效率
        W_comp0 = M1 * (H2 - H1)  # 压气机耗功
        

        #燃烧室计算        
        # dp_comb = 0.03
        # EFF_comb = 0.99
        # P3 = P2 * (1.0 - dp_comb)
        #质量分数
        # c1_m=0.91886
        # c2_m=0.05576
        # c3_m=0.0074
        # N2_m=0.00292
        # co2_m=0.01505
        # #获取组分
        # fuel=CoolProp.AbstractState('REFPROP','methane&ethane&propane&nitrogen&co2')
        # fuel.set_mass_fractions([c1_m, c2_m, c3_m, N2_m,co2_m])
        # molefrac_fuel = fuel.get_mole_fractions()

        # H_fuel = PropsSI('H', 'P', PF, 'T', TF, 'REFPROP::methane[{}]&ethane[{}]&propane[{}]&BUTANE[{}]&ISOBUTAN[{}]&PENTANE[{}]&IPENTANE[{}]&HEXANE[{}]&nitrogen[{}]&co2[{}]'.format(*molefrac_fuel))

        # air = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&argon&co2')
        # # air.set_mass_fractions([0.2314, 0.7552, 0.0129, 0.0005])
        # air.set_mass_fractions([0.22941, 0.74845, 0.020399, 1-0.22941-0.74845-0.020399])
        # molefrac_air = air.get_mole_fractions()

        # H0_fuel = PropsSI('H', 'T', 288.15, 'P', 101325, 'REFPROP::methane[{}]&ethane[{}]&propane[{}]&BUTANE[{}]&ISOBUTAN[{}]&PENTANE[{}]&IPENTANE[{}]&HEXANE[{}]&nitrogen[{}]&co2[{}]'.format(*molefrac_fuel))
        # H0_air = PropsSI('H', 'T', 288.15, 'P', 101325, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&argon[{}]&water[{}]'.format(*molefrac_air))
        # # frac = GasProp(patm, tdb, RH, M1, MF, c1, c2, c3, c4, ic4, c5, ic5, c6s, n2, o2, co2)
        # H0_gas = PropsSI('H', 'T', 288.15, 'P', 101325, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas))

        #H2=PropsSI('H','T',1257.6+273.15,'P',15.138*100000,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*frac))
        #LHV_fuel=((H2-H0_gas)*((SV1.MF + SVF.MF))-SV1.MF* (SV1.H-H0_air)-(Burn.H_fuel-H0_fuel)*SVF.MF)/SVF.MF
        # LHV_fuel=52197605.4 # ???

        #SV2.H = (SV1.MF* SV1.H+ SVF.MF*Burn.LHV_fuel+ Burn.H_fuel*SVF.MF) * Burn.EF/(SV1.MF + SVF.MF) 
        # H3 = (M1* (H2-H0_air)+ (H_fuel-H0_fuel)*MF+MF*LHV_fuel) * EFF_comb/(M1 + MF) +H0_gas
        #CB_P.EF = (SV2.H*SV2.MF)/(SV1.H*SV1.MF+SVF.MF*CB_P.LHV_fuel+CB_P.H_fuel*SVF.MF)
        
        # T3 = PropsSI('T','P',P3,'H',H3,'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas))

        #涡轮计算

        # S3 = PropsSI('S', 'T', T3, 'P', P3, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas))
        # H4_is = PropsSI('H', 'P', P4, 'S', S3, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas))
        # H4 = PropsSI('H', 'P', P4, 'T', T4, 'REFPROP::oxygen[{}]&nitrogen[{}]&co2[{}]&water[{}]&argon[{}]'.format(*molefrac_gas))

        # EFF_turb = (H3 - H4) / (H3 - H4_is)
        # W_turb = (M1 + MF) * (H3 - H4)

        W_GT0 = W_turb - W_comp0
        EFF_GT0 = W_GT0 / (MF * LHV_fuel)

        dataFrame0.loc[i-2] = [EFF_GT, EFF_GT0, W_GT, W_GT0]
        dataFrame0.to_csv('压损影响2 ' + localtime + '.csv')

        # print(i)

    return [EFF_comp, W_comp, EFF_turb, W_turb, EFF_GT, W_GT]

# 运行接口： 
result = run('C:/Users/ian_s/Desktop/Test/xnjs/input3.xlsx')