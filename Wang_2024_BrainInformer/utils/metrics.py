import numpy as np
from sklearn import metrics as me

def RSE(pred, true):
    # return np.sqrt(np.sum((true-pred)**2)) / np.sqrt(np.sum((true-true.mean())**2))
    return np.sqrt(((true-pred)**2).sum(0)) / np.sqrt(((true-true.mean())**2).sum(0))

def CORR(pred, true):
    u = ((true-true.mean(0))*(pred-pred.mean(0))).sum(0) 
    # d = np.sqrt(((true-true.mean(0))**2*(pred-pred.mean(0))**2).sum(0))
    d = np.sqrt(((true-true.mean(0))**2).sum(0))*np.sqrt(((pred-pred.mean(0))**2).sum(0))
    # return (u/d).mean(-1)
    return u/d

def CORR_sub(pred, true):
    u = ((true-true.mean(0))*(pred-pred.mean(0))).sum(0) 
    # d = np.sqrt(((true-true.mean(0))**2*(pred-pred.mean(0))**2).sum(0))
    d = np.sqrt(((true-true.mean(0))**2).sum(0))*np.sqrt(((pred-pred.mean(0))**2).sum(0))
    # return (u/d).mean(-1)
    return u/d
    # np.corrcoef(pred, true)
    # return 

def MAE(pred, true):
    # return np.mean(np.abs(pred-true))
    return (np.abs(pred-true)).mean(0)

def MSE(pred, true):
    # return np.mean((pred-true)**2)
    return ((pred-true)**2).mean(0)

def RMSE(pred, true):
    return np.sqrt(MSE(pred, true))
    # return np.sqrt(MSE(pred, true))

def MAPE(pred, true):
    return np.mean(np.abs((pred - true) / true))

def MSPE(pred, true):
    return np.mean(np.square((pred - true) / true))

def ACC(pred, true):
    return me.accuracy_score(true, pred)
def PRE(pred, true,labels):
    return me.precision_score(true, pred, labels,average='macro')
def RECALL(pred, true,labels):
    return me.recall_score(true, pred, labels, average='macro')
def F1(pred, true,labels):
    return me.f1_score(true, pred, labels, average='macro')
def CONFUSION(pred, true,labels):
    return me.confusion_matrix(true, pred, labels, average='macro')
def ROC(pred, true,labels):
    return me.roc_auc_score(true, pred)
def CLASS_REPORT(pred, true,labels,labelss):
    return me.classification_report(true, pred, labels, target_names=labelss)

def CONFUSION_MATRIX(pred, true):
    return me.confusion_matrix(true, pred)

def metric(pred, true):
    mae = MAE(pred, true)
    mse = MSE(pred, true)
    rmse = RMSE(pred, true)
    mape = MAPE(pred, true)
    mspe = MSPE(pred, true)
    rse = RSE(pred,true)
    corr = CORR(pred,true)
    
    return mae,mse,rmse,mape,mspe,rse,corr