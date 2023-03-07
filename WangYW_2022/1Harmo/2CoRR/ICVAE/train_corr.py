# -*- coding:utf-8 -*-
import vae_model as VAE
from vae_model import vae_model
import pandas as pd
import numpy as np
import scipy.io
import scipy.io
import os
from utils import onehot_test_label
import multiprocessing 
from math import ceil
subinfo = scipy.io.loadmat('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/SubInfo/SubInfo_420.mat')
site = subinfo['Site'].flatten()
print(type(site))
print(site.shape)

g_path = '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR'
workdir = ['Results','S2_Results']
metric_name = ['ALFF_FunImgARCW','fALFF_FunImgARCW','ReHo_FunImgARCWF','DegreeCentrality_FunImgARCWF','FC_D142']
#metric_name = ['fALFF_FunImgARCW_raw','DegreeCentrality_FunImgARCWF_raw']

def run():
    for wd in workdir:
        for name in metric_name:
        
            data_path = os.path.join(g_path,wd)
            print('reading %s %s ...' %(wd,name))
            features_struct = scipy.io.loadmat(os.path.join(data_path,"%s_raw.mat" % name))
            features = features_struct['raw']

            data = pd.DataFrame(features)
            

            labels = site
            Train = data
            Test = data
            train_label = labels
            test_label = np.full(site.shape[0],len(np.unique(site))-1)
            
            train_label = onehot_test_label(len(np.unique(site)),train_label)
            test_label = onehot_test_label(len(np.unique(site)),test_label)

            print(data.shape[1])
            
            pool=multiprocessing.Pool(5) 
            for i in range(ceil(data.shape[1]/512)): 
                path = os.path.join('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/VAE_related/new_vae_models',wd,name)
                if not os.path.exists(path):
                    os.makedirs(path) 
                    print('complete mkdir')
                h5_filename = os.path.join(path,"icvae_%i.h5"  % i)
                advh5 = os.path.join(path,"adv_%i.h5"  % i)
                
                if i != ceil(data.shape[1]/512)-1: 
                    train_data = Train.loc[:,i*512:i*512+511]
                    test_data = Test.loc[:,i*512:i*512+511]  
                else:
                    train_data = Train.iloc[:,-512:]
                    test_data = Test.iloc[:,-512:]

                pool.apply_async(vae_model,args=(train_data,train_label,10000,h5_filename,advh5,'train'))
                #vae_model(train_data,train_label,10000,h5_filename,advh5,'train')
            pool.close()
            pool.join()
            print(wd,name,'over') 

    

if __name__ == "__main__":
    #run()
    
    for wd in workdir:
        for name in metric_name:
            data_path = os.path.join(g_path,wd)
            print('reading %s %s ...' %(wd,name))
            features_struct = scipy.io.loadmat(os.path.join(data_path,"%s_raw.mat" % name))
            features = features_struct['raw']

            data = pd.DataFrame(features)
            print(data.shape)

            labels = site
            Train = data
            Test = data
            train_label = labels
            test_label = np.full(site.shape[0],len(np.unique(site))-1)
            print(test_label)
            train_label = onehot_test_label(len(np.unique(site)),train_label)
            test_label = onehot_test_label(len(np.unique(site)),test_label)
            all = pd.DataFrame()
            for i in range(ceil(data.shape[1]/512)):
                path = os.path.join('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/CoRR/VAE_related/new_vae_models',wd,name)
                if not os.path.exists(path):
                    os.makedirs(path)   
                
                h5_filename = os.path.join(path,"icvae_%i.h5"  % i)
                advh5 = os.path.join(path,"adv_%i.h5"  % i)
                
                #predict and save
                #map_label = onehot_test_label(3,np.full(data.shape[0],1))
                if i != ceil(data.shape[1]/512)-1:
                    p_data = data.loc[:,i*512:i*512+511]
                    x_hat = pd.DataFrame(VAE.vae_model(p_data,test_label,1,h5_filename,advh5,state='predict'))
                else:
                    p_data = data.iloc[:,-512:]
                    x_hat = pd.DataFrame(VAE.vae_model(p_data,test_label,1,h5_filename,advh5,state='predict'))
                    
                    cut = 512*i-data.shape[1] 
                    x_hat = x_hat.iloc[:,cut:]        
                    #print(x_hat.shape)
                all= pd.concat([all,x_hat],axis=1)
                 
            savepath = os.path.join('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/CORR/%s' % wd)
            if not os.path.exists(savepath):
                os.mkdir(savepath) 
            scipy.io.savemat(os.path.join(savepath,'%s_vae.mat' % name),{'ICVAE':all.to_numpy()})