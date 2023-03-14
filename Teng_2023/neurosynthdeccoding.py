#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 16:18:20 2022

@author: lihuixian
"""

import os
import nibabel as nib
from pprint import pprint
from nimare.extract import fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.dataset import Dataset
from nimare.decode import discrete

import numpy as np
from nilearn.plotting import plot_roi
from nimare.utils import get_resource_path

out_dir = os.path.abspath("./results_neurosynth_searchlight")
os.makedirs(out_dir, exist_ok=True)

dset = Dataset.load("./neurosynth_dataset.pkl.gz")

mask_filename = './Reslice_resultMask/3TwoExperimentOverlapresultNew.nii'
img = nib.load(mask_filename).get_fdata()

arr = np.zeros(dset.masker.mask_img.shape, int)
arr = arr + np.array(img, int)
mask_img = nib.Nifti1Image(arr, dset.masker.mask_img.affine)
# .. _neurosynth-roi-decoder-example:
#
# Decode an ROI image using the Neurosynth ROI association method
# -----------------------------------------------------------------------------

# This method decodes the ROI image directly, rather than comparing subsets of the Dataset like the
# other two.
decoder = discrete.ROIAssociationDecoder(mask_img)
decoder.fit(dset)

# The `transform` method doesn't take any parameters.
decoded_df = decoder.transform()

df = decoded_df.sort_values(by="r", ascending=False)
df.to_csv(os.path.join(out_dir, "3overlap_NeurosynthDecodingSearchlight.csv"))













