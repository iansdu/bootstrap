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

n_step = 120
n_pre = 24
train_pro =True
batch_size=13       #10-30
n_epochs = 100
learn_rate = 0.001
maxloops = 30
datalimit = [0.1,0.05,0.1,0.05]
#####train 
trainset=read_csv('traindata.csv',header=None)
trainset.to_csv('data.csv',index =None,header =None)
traindata_all=trainset.values
n_points = traindata_all.shape[0]
n_item = traindata_all.shape[1]
prediction = np.zeros((n_pre,n_item))
maxmin = np.zeros((n_item,2))
if train_pro:
    for prei in range(n_item):
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
 
        model=Sequential()

        #### 
        model.add(LSTM(512,input_shape=(train_x.shape[1],train_x.shape[2]),dropout=0))
        model.add(Dense(256,activation = 'sigmoid'))
        model.add(Dense(128,activation = 'sigmoid'))
        model.add(Dense(n_pre,activation = 'sigmoid'))
        adam = optimizers.adam(lr=learn_rate)
        model.compile(loss='mean_squared_error',optimizer=adam)
        model.fit(train_x,train_y,epochs=n_epochs,batch_size=batch_size,shuffle=True,verbose=1) 
        # verbose = 1 show cal process #shuffle - True disrupt the batch order
        model.save(os.path.join("Model" + str(prei+1)+ ".h5"))
        np.savetxt('maxmin.csv',maxmin, delimiter=",", fmt = '%s')