import numpy as np
import pandas as pd
from adv_model import train_adv
from utils import onehot_test_label
def adv_loss(raw_data,gen_data,advh5):
    y_true = np.full(raw_data.shape[0],1)
    y_true = onehot_test_label(3,y_true)
    y_false = np.full(gen_data.shape[0],0)
    y_false = onehot_test_label(3,y_false)
    
    fake_true = train_adv(raw_data,y_true,5,advh5,dis_trainable='False')
    fake_false = train_adv(gen_data,y_false,5,advh5,dis_trainable='False')
    adv_loss = np.mean(np.mean(fake_true)+ np.mean(fake_false))
    print('adv loss finished')
    return adv_loss
