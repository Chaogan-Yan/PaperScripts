#test
import vae_model as VAE
from vae_model import vae_model
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import scipy.io
import os
from utils import onehot_test_label
import multiprocessing 
from math import ceil

metric_name = ['ALFF_FunImgARCW','fALFF_FunImgARCW','ReHo_FunImgARCWF','DegreeCentrality_FunImgARCWF','FC_D142']

test_label = np.concatenate((np.full(6,1),np.full(6,1),np.full(6,1)),axis=0)
train_label = np.concatenate((np.full(41,0),np.full(41,1),np.full(41,2)),axis=0)

train_label = onehot_test_label(3,train_label)
test_label = onehot_test_label(3,test_label)


import multiprocessing
def run():  
    for name in metric_name:
        temp = name +'_raw.mat'
        features = scipy.io.loadmat('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS/'+temp)
        data = pd.DataFrame(features['alldata'])

        Train = pd.concat((data[0:41],data[41:82],data[82:123]),axis=0)
        Test = pd.concat((data[35:41],data[76:82],data[117:123]),axis=0)
        
        pool=multiprocessing.Pool(5) 
        for i in range(ceil(data.shape[1]/512)):
            path =os.path.join( '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSP3_ICVAE_ModelsParams/', name)
            if not os.path.exists(path):
                os.makedirs(path) 
                print('complete mkdir')
            h5_filename = os.path.join(path,"icvae_%i.h5"  % i)
            advh5 = os.path.join(path,"adv_%i.h5"  % i)
            
            # batch?
            if i != ceil(data.shape[1]/512)-1: 
                train_data = Train.loc[:,i*512:i*512+511]
                test_data = Test.loc[:,i*512:i*512+511]  
            else:
                train_data = Train.iloc[:,-512:]
                test_data = Test.iloc[:,-512:] 

            pool.apply_async(vae_model,args=(train_data,train_label,10000,h5_filename,advh5,'train'))
        pool.close()   
        pool.join()    
        


if __name__ == "__main__":
    #run()
 
    for name in metric_name:
        temp = name +'_raw.mat'
        features = scipy.io.loadmat('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS/'+temp)
        data = pd.DataFrame(features['alldata'])

        path = os.path.join( '/mnt/Data3/RfMRILab/Wangyw/harmonization_project/TST/TSP3_ICVAE_ModelsParams/', name)
        test_label = onehot_test_label(3,np.full(123,1)) 
        #print(test_label)
        data_hat = pd.DataFrame()
        # find the corresponding model
        for i in range(ceil(data.shape[1]/512)):
            icvae_h5 = os.path.join(path,"icvae_%i.h5"  % i)
            adv_h5 = os.path.join(path,"adv_%i.h5"  % i)
            #predict and save         
            if i != ceil(data.shape[1]/512)-1:
                p_data = data.loc[:,i*512:i*512+511]
                x_hat = pd.DataFrame(VAE.vae_model(p_data,test_label,1,icvae_h5,adv_h5,state='predict'))
            else:
                p_data = data.iloc[:,-512:]
                x_hat = pd.DataFrame(VAE.vae_model(p_data,test_label,1,icvae_h5,adv_h5,state='predict'))
                
                cut = 512*i-data.shape[1] 
                x_hat = x_hat.iloc[:,cut:]        
                

            data_hat= pd.concat([data_hat,x_hat],axis=1)
        print(data_hat)    
        temp = name +'_ICVAE.mat'
        scipy.io.savemat(os.path.join('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS',temp),{'ICVAE':data_hat.to_numpy()})
        # data_hat.to_csv(os.path.join('/mnt/Data3/RfMRILab/Wangyw/harmonization_project/Restart/HarmonizationResults/TSP3/ResultsS',temp),header=False,index=False)
