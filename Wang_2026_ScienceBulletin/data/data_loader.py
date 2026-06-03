from audioop import cross
from ctypes import sizeof
import os
import sys
import numpy as np
import pandas as pd
import h5py  
from numpy.core.shape_base import hstack, vstack
import random
from utils.checkdata import check_data_nan, check_data_zero
from silx.io.dictdump import h5todict
import time
from sklearn.model_selection import KFold
from utils.read_excel import get_label
from utils.log_utils import get_logger
from utils.func import pearson
LOGGER = get_logger("informer")

import torch
from torch.utils.data import Dataset, DataLoader

from utils.tools import StandardScaler
from utils.timefeatures import time_features
from utils.rolling import rolling_window

import warnings
warnings.filterwarnings('ignore')

class Dataset_Multi_Task(Dataset):
    def __init__(self, root_path, flag='train', size=None, 
                 features='S', data_path='UKBiobank_Group1.csv', 
                 target='OT', scale=True, inverse=False, timeenc=0, freq='h', cols=None, data_fine_tune = '', setting = '', cross = cross): 
        if size == None:
            self.seq_len = 24*4*4
            self.label_len = 24*4
            self.pred_len = 24*4
        else:
            self.seq_len = size[0]
            self.label_len = size[1]
            self.pred_len = size[2]
        # init
        assert flag in ['train', 'test', 'val', 'all']
        type_map = {'train':0, 'val':1, 'test':2, 'all':3}
        self.set_type = type_map[flag]
        
        self.features = features
        self.target = target
        self.scale = scale
        self.inverse = inverse
        self.timeenc = timeenc
        self.freq = freq
        self.cols=cols
        self.root_path = root_path
        self.data_path = data_path
        self.len_x=list()
        self.len_y=list()
        self.num_win = list()
        self.len_data = list()
        self.data_num = 0
        self.data_index = list()
        self.cross = cross 
        self.data_fine_tune = data_fine_tune
        self.setting = setting
        self.data_label = list()
        self.data_station = list()
     

        self.__read_data__()

    def __read_data__(self):
        self.scaler = StandardScaler()
        f = h5py.File(os.path.join(self.root_path,
                                        self.data_path),'r') 
        if self.data_path == "MDD_4task.hdf5":
            label_file_name = 'MDD_4task.xls'
        elif self.data_path == "MULTI_RUMI_EMO.hdf5":
            label_file_name = 'rumi_emo.xls'
        elif self.data_path == "Rumination_Suzhou.hdf5":
            label_file_name = 'Rumination_Suzhou.xls'
        task_label, station_label = get_label(label_file_name,self.root_path,0,1,0,0)
        data_num = int(len(f))
        num_train = int(data_num*0.9)
        num_test = int(data_num*0.1)
        num_vali = int(data_num - num_train - num_test) 

        index_random = np.arange(data_num)
        random.shuffle(index_random)

        kf_tt = KFold(n_splits=10,shuffle=False) 
        kf_tv = KFold(n_splits=9,shuffle=False)
        train_cross = []
        test_cross = []
        valid_cross = []
        split1 = 0

        for trainvalid_index, test_index in kf_tt.split(index_random):
            split2 = 0
            for train_index , valid_index in kf_tv.split(index_random[trainvalid_index]):
                if split1 == split2 or split1==9: 
                    train_cross.append(index_random[trainvalid_index][train_index].tolist() ) 
                    valid_cross.append(index_random[trainvalid_index][valid_index].tolist() )
                    test_cross.append(index_random[test_index].tolist() )
                    break
                split2 = split2 + 1
            split1 = split1 + 1
        cross_index = {0:train_cross[self.cross], 1:valid_cross[self.cross], 2:test_cross[self.cross], 3:np.arange(data_num)}
        data_index = cross_index[self.set_type]

        data_name=list()
        for key in f.keys():
            data_name.append(key)

        k = 0
        for di in data_index:
            keys = data_name[di]
            dataforlen = f[keys][()]
            self.len_data.append(len(dataforlen))
            self.num_win.append(len(dataforlen)-self.seq_len-self.pred_len+1)
            self.len_x.append(self.num_win[k]*self.seq_len)
            self.len_y.append(self.num_win[k]*(self.label_len+self.pred_len))
            k = k+1
        self.data_num = int(k)
        if self.features=='M' or self.features=='MS':
            col = len(dataforlen[0])
        elif self.features=='S':
            col = 1
        dataMix_x = np.zeros((sum(self.len_x),col),dtype = None, order = 'C')
        dataMix_y = np.zeros((sum(self.len_y),col),dtype = None, order = 'C')
        data_label = np.zeros((sum(self.len_x)),dtype = None, order = 'C')
        data_station = np.zeros((sum(self.len_x)),dtype = None, order = 'C')
   
        i = 0
        self.data_index = data_index
        for di in data_index:
            keys = data_name[di]
        
            t_before = time.time()
            data_single = f[keys][()] 
            data_label_sub = task_label[di]
            data_station_sub = station_label[di]

            time_point = len(data_single)
            data_window_num = time_point-self.seq_len- self.pred_len + 1 

            if self.features=='M' or self.features=='MS':
                cols_num= len(data_single[0])
                
            
            elif self.features=='S':
                data_single = data_single[:,self.target]
                cols_num= 1

            if self.scale:
                train_data = data_single
                self.scaler.fit(data_single)
                data = self.scaler.transform(data_single) 
            else:
                data = data_single
            data_single_x = data 

            if self.inverse: 
                data_single_y = data_single
            else:
                data_single_y = data
       


            data_split_x = np.reshape(data_single_x,(time_point,cols_num)) 
            data_split_y = np.reshape(data_single_y,(time_point,cols_num))
            try:
                seq_x = rolling_window(data_split_x[:data_split_x.shape[0]-self.pred_len],self.seq_len, self.pred_len)
                seq_y = rolling_window(data_split_y[(self.seq_len-self.label_len):],(self.label_len+self.pred_len), self.pred_len)
            except ValueError:
                print(keys, di, data_split_x.shape, data_split_y.shape)
            seq_x = seq_x.reshape(seq_x.shape[0]*seq_x.shape[1],seq_x.shape[2])
            seq_y = seq_y.reshape(seq_y.shape[0]*seq_y.shape[1],seq_y.shape[2])
            
            if np.isnan(seq_x).any() == 1:
                check_x_single = np.argwhere(np.isnan(seq_x))
            else:
                check_x_single = 1
            try:
                if check_x_single != 1:
                    print(keys)
                    print(check_x_single)
            except ValueError:
                print(keys)

            dataMix_x[sum((self.len_x[0:i])):sum(self.len_x[0:i+1]),] = seq_x
            dataMix_y[sum(self.len_y[0:i]):sum(self.len_y[0:i+1]),] = seq_y
            data_label[sum((self.len_x[0:i])):sum(self.len_x[0:i+1])] = data_label_sub
            data_station[sum((self.len_x[0:i])):sum(self.len_x[0:i+1])] = data_station_sub


         
            i = i+1
          

        self.seq_x = dataMix_x
        self.seq_y = dataMix_y
        self.seq_label = data_label.astype(np.int64)
        self.seq_station = data_station.astype(np.int64)
        LOGGER.info("data_num:{}".format(self.data_num))




        f.close()  

    def __getitem__(self, index):
        seq_x = self.seq_x[index*self.seq_len:(index+1)*self.seq_len]
        seq_y = self.seq_y[index*(self.pred_len+self.label_len):(index+1)*(self.pred_len+self.label_len)]
        seq_label = self.seq_label[index*self.seq_len:(index+1)*self.seq_len][0]
        seq_station = self.seq_station[index*self.seq_len:(index+1)*self.seq_len][0]
        return seq_x, seq_y, seq_label, seq_station
    
    def __len__(self):
        return int(len(self.seq_x)/self.seq_len) 


    def inverse_transform(self, data):
        return self.scaler.inverse_transform(data)