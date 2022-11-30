#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 19:24:53 2021

@author: Lihuixian
"""

import os
import glob
import numpy as np
import pandas as pd
import pingouin as pg
import nibabel as nib
from neurora.stuff import get_affine  
from neurora.rdm_cal import fmriRDM
from neurora.corr_cal_by_rdm import fmrirdms_corr
from neurora.nii_save import corr_save_nii

workingpath = "/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/MeanBlodAnalysis/RSAnewprepocess"
fmriRDMsOutpath = os.path.join(workingpath, "section1", "fmriRDMresultSearchlight5")
RSAresultpath = os.path.join(workingpath, "section1", "partialSpearman_RSAresults5")
RSAresultZpath = os.path.join(workingpath, "section1", "partialSpearman_RSAresultsZ5")
bhvpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RSAnewprocess/bhvRDMresult/csvfiles'
bhvsentenceLpath='/mnt/Data3/RfMRILab/Lihuixian/DataAnalysis/TaskAnalysis/2020WYSWYT/RSAnewprocess/partialSpearmanAnalysis/bhvRDMsentenceLength/csvfiles'

"""*** Calculating the RDM by Searchlight ***"""
mask_filename = "/mnt/Data3/RfMRILab/Lihuixian/DPABI_V6.0_ForCamp/Templates/BrainMask_05_61x73x61.img"
maskdata = nib.load(mask_filename).get_fdata()
nx, ny, nz = maskdata.shape
datapath = os.path.join(workingpath, "section1", "MeanBlodSpeakEvent")
subID = os.listdir(datapath)
for sub in range(len(subID)):
    subname = os.path.join(datapath, subID[sub])
    print(subname)
    speakevent = glob.glob("%s/*.nii" % subname)
    ncon = len(speakevent)
    fmri_data = np.full([ncon, nx, ny, nz], np.nan)
    for event in range(len(speakevent)):
        fmri_data[event] = nib.load(speakevent[event]).get_fdata()
    ## reshaoe the data: [ncon, nx, ny, nz] -> [ncon, nsubs, nx, ny, nz]
    ## here just one subject's data
    fmri_data = np.reshape(fmri_data, [ncon, 1, nx, ny, nz])
    print(fmri_data.shape)
    ## calculate the RDMs by Searchlight voxels=5 b=1 # here just one subject's data sub_opt=0
    print("calculate fmri_RDMs...")
    fmri_RDMs = fmriRDM(
        fmri_data,
        ksize=[5, 5, 5],
        strides=[1, 1, 1],
        sub_opt=0,
        method="correlation",
        abs=True,
    )
    print(fmri_RDMs.shape)
    ## save fmri_RDMs
    savepath = "%s/fmriRDMs_%s.npy" % (fmriRDMsOutpath, subID[sub])
    np.save(savepath, fmri_RDMs)
    print(subID[sub] + "The fmri_RDMs was saved successfully!")

    """***Calculating the representational similarities***"""
    bhvsubpath = "%s/RDM_%s.csv" % (bhvpath, subID[sub])
    bhv_RDM = np.loadtxt(bhvsubpath, delimiter=",")
    # covariation：based on sentence length
    bhvsentenceLpathsub = os.path.join(
        bhvsentenceLpath, "wordNumRDM_%s.csv" % (subID[sub])
    )
    bhvsentenceL_RDM = np.loadtxt(bhvsentenceLpathsub, delimiter=",")
    print("calculate partial spearman RSA..." + subID[sub])
    # calculate the number of the calculation units in the x, y, z directions
    n_x = np.shape(fmri_RDMs)[0]
    n_y = np.shape(fmri_RDMs)[1]
    n_z = np.shape(fmri_RDMs)[2]

    # initialize the corrs
    corrs_result = np.full([n_x, n_y, n_z, 2], np.nan)
    Zcorrs_result = np.full([n_x, n_y, n_z, 2], np.nan)
    # get number of conditions
    cons = np.shape(bhv_RDM)[0]
    n = int(cons * (cons - 1) / 2)
    print(np.shape(bhv_RDM))
    print(np.shape(fmri_RDMs))
    print(np.shape(bhvsentenceL_RDM))
    v1_bhvRDM = np.zeros([n], dtype=np.float64)
    v2_fmriRDM = np.zeros([n], dtype=np.float64)
    v3_bhvsentenceLRDM = np.zeros([n], dtype=np.float64)

    # assignment
    nn = 0
    for i in range(cons - 1):
        for j in range(cons - 1 - i):
            v1_bhvRDM[nn] = bhv_RDM[i, i + j + 1]
            v3_bhvsentenceLRDM[nn] = bhvsentenceL_RDM[i, i + j + 1]
            nn = nn + 1

    # calculate the corrs: partial spearman
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                fmri_rdms = fmri_RDMs[i, j, k]
                # get fmri voxel RDM vector
                mm = 0
                for ifmri in range(cons - 1):
                    for jfmri in range(cons - 1 - ifmri):
                        v2_fmriRDM[mm] = fmri_rdms[ifmri, ifmri + jfmri + 1]
                        mm = mm + 1
                # sort data to satisfy format
                data = pd.DataFrame(
                    {
                        "bhv_RDM": v1_bhvRDM,
                        "fmri_RDMs": v2_fmriRDM,
                        "Covariation": v3_bhvsentenceLRDM,
                    }
                )
                # nan-->0
                dataxin = data.fillna(0)
                corrs = pg.partial_corr(
                    dataxin,
                    x="bhv_RDM",
                    y="fmri_RDMs",
                    covar="Covariation",
                    method="spearman",
                )
                corrs_result[i, j, k] = corrs[["r", "p-val"]]

                ##Do the Fisher-Z transform of the r
                Zcorrs_result[i, j, k] = corrs[["r", "p-val"]]
                Zcorrs_result[i, j, k, 0] = 0.5 * np.log(
                    (1 + corrs["r"]) / (1 - corrs["r"])
                )
    """*** sav RSA result ***"""
    # get the affine info
    affine = get_affine(mask_filename)
    # save the RSA result as a .nii file
    RSAresultfilename = "%s/partialSpearmanRSAimg_%s.nii" % (
        RSAresultpath,
        subID[sub],
    )
    # If img_background=None, the background will be ch2.nii.gz.
    img = corr_save_nii(
        corrs_result,
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

    RSAZresultfilename = "%s/ZpartialSpearmanRSAimg_%s.nii" % (
        RSAresultZpath,
        subID[sub],
    )
    imgZ = corr_save_nii(
        Zcorrs_result,
        filename=RSAZresultfilename,
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
