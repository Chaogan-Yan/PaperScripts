
import numpy as np
def pearson(x):
    pearson_corr_x = (np.corrcoef(x,rowvar=False))
    triu_x = pearson_corr_x[np.triu_indices(pearson_corr_x.shape[0], k=1)] 
    return triu_x

def pearson_only(x):
    pearson_corr_x = (np.corrcoef(x,rowvar=False))
    return pearson_corr_x

def softmax(x):
    """ softmax function """
    
    # assert(len(x.shape) > 1, "dimension must be larger than 1")
    # print(np.max(x, axis = 1, keepdims = True)) # axis = 1, 行
    
    x -= np.max(x, axis = 1, keepdims = True) 
        
    x = np.exp(x) / np.sum(np.exp(x), axis = 1, keepdims = True)
    
    return x