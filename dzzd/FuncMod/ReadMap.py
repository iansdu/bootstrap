# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:29 2019

@author: dell
"""
import numpy
from scipy import interpolate
import math
import scipy

def ReadMap(CN_T, PR_T,SZ, map,CM_T = 0):
    if CN_T > map.AX[0] and CN_T <map.AX[-1]:
        C_flag=0 #实际转速在map里
    else:
        C_flag=1 #实际转速超出map
        print(C_flag)
        if CN_T < map.AX[0]:
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