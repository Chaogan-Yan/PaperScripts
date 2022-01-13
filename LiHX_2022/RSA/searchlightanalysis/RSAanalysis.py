#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 19:24:53 2021

@author: Lihuixian
"""

import os
import glob
import numpy as np
import nibabel as nib
from neurora.stuff import get_affine #datamask ROIanalysis
#from neurora.rsa_plot import plot_rdm
from neurora.rdm_cal import  fmriRDM
from neurora.corr_cal_by_rdm import fmrirdms_corr
from neurora.nii_save import corr_save_nii

workingpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess'
fmriRDMsOutpath=os.path.join(workingpath,'section2','fmriRDMresultSearchlight5')
RSAresultpath=os.path.join(workingpath,'section2','RSAresults5')
RSAresultZpath=os.path.join(workingpath,'section2','RSAresultsZ5')
bhvpath=os.path.join(workingpath,'bhvRDMresult','csvfiles')
"""*** Calculating the RDM by Searchlight ***"""
mask_filename='/mnt/Data3/RfMRILab/Lihuixian/DPABI_V5.1_201230/Templates/BrainMask_05_61x73x61.img'
maskdata = nib.load(mask_filename).get_fdata()
nx, ny, nz = maskdata.shape
datapath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/WYSWYT/RSAanalysisWYSWYT/RSAnewprepocess/section2/SpeakEventBeta'
subID=os.listdir(datapath)
for sub in range(len(subID)):
    subname=os.path.join(datapath,subID[sub])
    print(subname)
    speakeventbeta=glob.glob('%s/*.nii' %subname)
    ncon=len(speakeventbeta)
    fmri_data = np.full([ncon, nx, ny, nz], np.nan)
    for event in range(len(speakeventbeta)):
        fmri_data[event] = nib.load(speakeventbeta[event]).get_fdata()
    fmri_data = np.reshape(fmri_data, [ncon,1, nx, ny, nz])
    print(fmri_data.shape)
    ## calculate the RDMs by Searchlight voxels=5 b=1
    print('calculate fmri_RDMs...')
    fmri_RDMs = fmriRDM(fmri_data,ksize=[5, 5, 5],strides=[1, 1, 1],sub_opt=0,method='correlation', abs=True)
    print(fmri_RDMs.shape)
    ## save fmri_RDMs
    savepath=('%s/fmriRDMs_%s.npy' %(fmriRDMsOutpath,subID[sub]))
    np.save(savepath,fmri_RDMs)
    print(subID[sub]+'The fmri_RDMs was saved successfully!')
    """***Calculating the representational similarities***"""
    bhvsubpath=('%s/RDM_%s.csv' %(bhvpath,subID[sub]))
    bhv_RDM = np.loadtxt(bhvsubpath,delimiter=",")
    print('calculate RSA...'+subID[sub])         
    corrs=fmrirdms_corr(bhv_RDM,fmri_RDMs,method="spearman", fisherz=False, rescale=False, permutation=False, iter=5000)
    print('calculate ZRSA...'+subID[sub]) 
    #fisherz_rdm;not permutation test
    corrsZ=fmrirdms_corr(bhv_RDM,fmri_RDMs,method="spearman", fisherz=True, rescale=False, permutation=False, iter=5000)
    """*** sav RSA result ***"""
    # get the affine info
    affine = get_affine(mask_filename)
    # save the RSA result as a .nii file
    RSAresultfilename = ('%s/RSAimg_%s.nii' %(RSAresultpath,subID[sub]))
    #If img_background=None, the background will be ch2.nii.gz.
    img = corr_save_nii(corrs, filename=RSAresultfilename,affine=affine,corr_mask=mask_filename,size=[61, 73, 61],ksize=[5, 5, 5],strides=[1, 1, 1], p=1,r=0,correct_method=None,smooth=False,plotrlt=False,img_background=None)
    
    RSAZresultfilename = ('%s/ZRSAimg_%s.nii' %(RSAresultZpath,subID[sub]))
    imgZ = corr_save_nii(corrsZ, filename=RSAZresultfilename,affine=affine,corr_mask=mask_filename,size=[61, 73, 61],ksize=[5, 5, 5],strides=[1, 1, 1], p=1,r=0,correct_method=None,smooth=False,plotrlt=False,img_background=None)





