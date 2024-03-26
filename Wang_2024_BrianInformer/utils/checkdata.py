
import sys
sys.path.append("/mnt/Data5/RfMRILab/Wangzh/Code/informer/Informer2020-main/")
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from utils.metrics import metric,CORR
import pandas as pd
import h5py
import os
 
 
def check_data_nan(data_check):
    if np.isnan(data_check).any() == 1:
        data_nan = np.argwhere(np.isnan(data_check))
        return data_nan
    return 1
def check_data_zero(data_check):
    if np.any(data_check == 0):
        data_zero = np.argwhere(np.any(data_check == 0))
        return data_zero
    return 1