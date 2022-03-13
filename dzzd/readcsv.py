# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 16:55:41 2019

@author: dell
"""
import csv
import numpy
import xlwt

#  将数据写入新文件
def data_write(file_path, datas):
    f = xlwt.Workbook()
    sheet1 = f.add_sheet(u'sheet1',cell_overwrite_ok=True) #创建sheet
    
    #将数据写入第 i 行，第 j 列
    i = 0
    for data in datas:
        for j in range(len(data)):
            sheet1.write(i,j,data[j])
        i = i + 1
        
    f.save(file_path) #保存文件

def transpose(self, matrix):
        new_matrix = []
        for i in range(len(matrix[0])):
            matrix1 = []
            for j in range(len(matrix)):
                matrix1.append(matrix[j][i])
            new_matrix.append(matrix1)
        return new_matrix

file = "C:/Users/dell/Desktop/西部管道功率测算/霍尔果斯气路1112/h4/h4.csv"
with open(file,'r', encoding = 'UTF-8') as f :
    # 使用csv.DictReader读取文件中的信息
    reader = csv.DictReader(f)
    N1 = []
    N2 = []
    T1 = []
    P1 = []
    T2 = []
    P2 = []
    T34 = []
    P34= []
    T4 = []
    P4 = []
    P2f = []
    P1f = []
    Tf = []
    Pf = []
    id_1 = []
    time_1 = []
    Ctin = []
    Cpin = []
    Ctou = []
    Cpou = []
    for row in reader :
        if (row['time'] == '2019/11/8 12:33' or 
            row['time'] == '2019/11/9 12:22' or
            row['time'] == '2019/11/9 15:03' or
            row['time'] == '2019/11/10 16:21' or
            row['time'] == '2019/11/11 11:31' or
            row['time'] == '2019/11/11 11:55' or
            row['time'] == '2019/11/11 15:13' or
            row['time'] == '2019/11/11 15:24'):
            time_1.append(str(row['time']))
            id_1.append(float(row['id']))
            N1.append(float(row['T1_NGGSEL']))
            N2.append(float(row['T1_NPTSEL']))
            T1.append(float(row['T1_T2SEL']))
            P1.append(float(row['T1_P2SEL']))
            T2.append(float(row['T1_T3SEL']))
            P2.append(float(row['T1_PS3SEL']))
            T34.append(float(row['T1_T48SEL']))
            P34.append(float(row['T1_HWIN_P48A']))
            T4.append(float(row['T1_T8SEL']))
            #P4.append(float(row['T1_T8SEL']))
            P4.append(float(101325))
            P2f.append(float(row['T1_gp2a']))
            P1f.append(float(row['T1_gp1a']))
            Tf.append(float(row['T1_tfuela']))
            Pf.append(float(row['T1_gfmvpsfbk']))
            Ctin.append(float(row['T1_a26gs']))
            Cpin.append(float(row['T1_a63gs']))
            Ctou.append(float(row['T1_a26gda']))
            Cpou.append(float(row['T1_a63gd']))
        else:
            pass
    datadiag = []  
    data = (id_1,time_1,N1,N2, T1, P1, T2, P2, T34, P34, T4, P4, P2f, P1f,Tf, Pf,Ctin,Cpin,Ctou,Cpou)
    data = transpose((),data)
    n = len(data)
    for i in range(n):
        datadiag.append(data[n-i-1])
    data_write('datadiag.xls',datadiag)
    #numpy.savetxt('datadiag.csv', datadiag,delimiter = ',')
#
        #转变数据内形式，str2int
        #home_team_goals.append(int(row['Home Team Goals']))
        #away_team_goals.append(int(row['Away Team Goals']))

 
