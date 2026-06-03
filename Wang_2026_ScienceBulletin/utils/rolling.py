import numpy as np
import pandas as pd
def rolling_window(a, window, pred_len):
    shape = np.array([a.shape[0]-window+1,window,a.shape[1]])
    strides = np.array([a.shape[1]*8,a.shape[1]*8,8]) 
    result = np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
    return result
def rolling1(a, window):
    shape = (a.size - window + 1, window)
    strides = (a.itemsize, a.itemsize)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
def rolling(a, window):
    df = pd.DataFrame(a)
    return df.rolling(window, min_periods=1)
