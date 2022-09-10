#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 19:24:53 2021

@author: Lihuixian
"""

import os
import numpy as np
from neurora.stuff import get_affine 
from neurora.corr_cal_by_rdm import fmrirdms_corr
from neurora.nii_save import corr_save_nii

workingpath = "/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/RandomAnalysis"
RSAresultpath = os.path.join(workingpath, "random_RSAresults5")
bhvpath = "/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RSAnewprocess/bhvRDMresultRandom/csvfiles"

mask_filename = "/mnt/Data3/RfMRILab/Lihuixian/DPABI_V6.0_ForCamp/Templates/BrainMask_05_61x73x61.img"
# have calculated nerual rdm
datapath = "/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/RSAnewprepocess/fmriRDMresultSearchlight5"

subpath = "/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/MeanBlodSpeakEvent47"
subID = os.listdir(subpath)
for sub in range(len(subID)):
    subname = os.path.join(datapath, "fmriRDMs_%s.npy" % (subID[sub]))
    print(subname)
    fmri_RDMs = np.load(subname)
    """***Calculating the representational similarities***"""
    bhvsubpath = os.path.join(bhvpath, "randomRDM_%s.csv" % (subID[sub]))
    bhv_RDM = np.loadtxt(bhvsubpath, delimiter=",")
    print("calculate RSA..." + subID[sub])
    corrs = fmrirdms_corr(
        bhv_RDM,
        fmri_RDMs,
        method="spearman",
        fisherz=False,
        rescale=False,
        permutation=False,
        iter=5000,
    )
    print("calculate ZRSA..." + subID[sub])
    """*** sav RSA result ***"""
    # get the affine info
    affine = get_affine(mask_filename)
    # save the RSA result as a .nii file
    RSAresultfilename = "%s/randomRSAimg_%s.nii" % (RSAresultpath, subID[sub])
    # If img_background=None, the background will be ch2.nii.gz.
    img = corr_save_nii(
        corrs,
        filename=RSAresultfilename,
        affine=affine,
        corr_mask=mask_filename,
        size=[61, 73, 61],
        ksize=[5, 5, 5],
        strides=[1, 1, 1],
        p=1,
        r=-1,
        correct_method=None,
        smooth=False,
        plotrlt=False,
        img_background=None,
    )
