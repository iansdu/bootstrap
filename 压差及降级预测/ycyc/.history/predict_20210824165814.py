import keras
from math import sqrt
from numpy import concatenate
from matplotlib import pyplot
import pandas as pd
from pandas import read_csv
from pandas import DataFrame
from pandas import concat
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import mean_squared_error
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.callbacks import TensorBoard
from IPython.display import SVG
from keras.utils.vis_utils import model_to_dot
from keras.layers import GRU, Dense, Masking, Dropout, Activation
from sklearn import preprocessing   
import numpy as np
from keras.models import load_model  
import h5py
from keras.models import Model #泛型模型
from keras.layers import Input
import os
from keras import optimizers
import matplotlib.pyplot as plt
from scipy import optimize

n_step = 120
n_pre = 24
train_pro =False
batch_size=13       #10-30
n_epochs = 100
learn_rate = 0.001
maxloops = 100
datalimit = [0.1,0.05,0.1,0.05]
#####train 
trainset=read_csv('prePredict.csv',header=None)
trainset.to_csv('data.csv',index =None,header =None)
traindata_all=trainset.values
n_points = traindata_all.shape[0]
n_item = traindata_all.shape[1]
prediction = np.zeros((n_pre,n_item))
maxmin = np.zeros((n_item,2))

model = load_model(os.path.join("Model"+str(0+1)+ ".h5"))

for nlp in range(maxloops):
#for prei in range(n_item):
    trainset=read_csv('data.csv',header=None)
    trainsetorg = read_csv('data.csv',header=None)
    traindata_all=trainset.values
    n_points = traindata_all.shape[0]
    n_item = traindata_all.shape[1]
    for prei in range(n_item):
    #for nlp in range(maxloops):    
        traindata = traindata_all[:,prei]
        traindata = traindata.reshape((-1,1))
        data_max = np.max(traindata)
        data_min = np.min(traindata)
        maxmin[prei,0] = data_max
        maxmin[prei,1] = data_min
        n_data = n_points-n_step-n_pre
        for n_in in range(n_points):
            traindata[n_in,:] = (traindata[n_in,:] - data_min) / (data_max-data_min)
    
        train_x = np.zeros((n_data,n_step+1))
        train_y = np.zeros((n_data,n_pre))
        for i in range(n_data):
            for k in range(n_pre):
                train_y[i,k] = traindata[i+n_step+k,0]
            for j in range(n_step):
                train_x[i,j]= traindata[i+j,0]
                train_x[i,n_step]=i/2000
    
        train_x=train_x.reshape((n_data,1,n_step+1))
        train_y=train_y.reshape((n_data,n_pre))
        initialization = 'glorot_normal'
 
        #model=Sequential()
        ####unnote this part
        #model = load_model(os.path.join("BanshanModelDEC"+str(prei+1)+ ".h5")) # replace training part with this sentence if you have trained a model
        maxmin=read_csv('maxmin.csv',header=None)
        maxmin=maxmin.values
        data_max = maxmin[prei,0]
        data_min = maxmin[prei,1]
        test_x = np.zeros((1,1,n_step+1))
        test_y = np.zeros((1,n_pre))
        pre_y = np.zeros((n_pre))
        for i in range(n_step):
            test_x[0,0,i] = traindata[n_points-n_step+i]
        
        test_y[0] = model.predict(test_x)
        for i in range(n_pre):
            pre_y[i] = test_y[0,i] *(data_max-data_min) + data_min
        prediction[:,prei] = pre_y
    prediction_df = pd.DataFrame(prediction)
    trainsetorg=pd.concat([trainsetorg,prediction_df])
    np.savetxt('prediction.csv', prediction, delimiter=",", fmt = '%s')
    trainsetorg.to_csv('data.csv',index =None,header =None)

    print(nlp)

def f1(x,a,b):
        return a*x+b#**2+b*x+c

DP=trainsetorg.values[:,0]*100
x=np.arange(1,DP.size+1)
# f1=np.polyfit(x,DEC,2)
# p1=np.poly1d(f1)
# DEC_trend=p1(DEC)
A1,B1=optimize.curve_fit(f1,x,DP)[0]
DP_trend=A1*x+B1#**2+B1*x+C1
#if DEC_trend[-1]-DEC_trend[0]>2.5 :  #效率降级大于2%
# if DP[-1]>2:
plot1=plt.plot(x,DP_trend,'r')
plot2=plt.plot(x,DP,'s')
plt.show()
    # break


# for nlp in range(maxloops):

#     trainset=read_csv('prePredict.csv',header=None)
#     trainsetorg = read_csv('data.csv',header=None)
#     traindata_all=trainset.values
#     n_points = traindata_all.shape[0]
#     n_item = traindata_all.shape[1]
#     for prei in range(n_item):
#         traindata = traindata_all[:,prei]
#         traindata = traindata.reshape((-1,1))
#         # data_max = np.max(traindata)
#         # data_min = np.min(traindata)
#         # maxmin[prei,0] = data_max
#         # maxmin[prei,1] = data_min
#         # n_data = n_points-n_step-n_pre
#         # for n_in in range(n_points):
#         #     traindata[n_in,:] = (traindata[n_in,:] - data_min) / (data_max-data_min)
    
#         # train_x = np.zeros((n_data,n_step))
#         # train_y = np.zeros((n_data,n_pre))
#         # for i in range(n_data):
#         #     for k in range(n_pre):
#         #         train_y[i,k] = traindata[i+n_step+k,0]
#         #     for j in range(n_step):
#         #         train_x[i,j]= traindata[i+j,0]
    
#         # train_x=train_x.reshape((n_data,1,n_step))
#         # train_y=train_y.reshape((n_data,n_pre))
#         # initialization = 'glorot_normal'
 
#         model=Sequential()
#         ####unnote this part
#         model = load_model(os.path.join("Model"+str(prei+1)+ ".h5")) # replace training part with this sentence if you have trained a model
#         maxmin=read_csv('maxmin.csv',header=None)
#         maxmin=maxmin.values
#         data_max = maxmin[prei,0]
#         data_min = maxmin[prei,1]
#         test_x = np.zeros((1,1,n_step))
#         test_y = np.zeros((1,n_pre))
#         pre_y = np.zeros((n_pre))
#         for i in range(n_step):
#             test_x[0,0,i] = traindata[n_points-n_step+i]
#         test_y[0] = model.predict(test_x)
#         for i in range(n_pre):
#             pre_y[i] = test_y[0,i] *(data_max-data_min) + data_min
#         prediction[:,prei] = pre_y
#     prediction_df = pd.DataFrame(prediction)
#     trainsetorg=pd.concat([trainsetorg,prediction_df])
#     np.savetxt('prediction.csv', prediction, delimiter=",", fmt = '%s')
#     #trainsetorg.to_csv('data.csv',index =None,header =None)

