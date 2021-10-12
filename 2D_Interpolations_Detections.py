#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:42:42 2021

@author: Sherry Wong
"""

#%%
## Windows Filenames
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/wongb/Documents/Python_Scripts"
#%%
## Mac filenames
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/bwong/URS_FLASH_DataParker"

#%%

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker, cm #for log scale
sys.path.insert(0, pwd)
# import copy

import athena_read as ath
import h5py
import hdf5_parser

import scipy.interpolate as interp

#%%
newOut = hdf5_parser.setup(filename, format = "cartesian")
x = newOut['posXarray'].flatten()
y = newOut['posYarray'].flatten()
z = newOut['densityArray'].flatten()
# from matplotlib import ticker, cm #for log scale
# plt.tricontourf(x, y, z, levels=100, locator=ticker.LogLocator())

lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
plt.clf()
plt.tricontourf(x, y, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
plt.savefig(pwd + "/Plots/ZcleanB1.png")

#%% interpolate

x = newOut['magXarray'][0]
y = newOut['posYarray'][0]
Xslice = interp.interp1d(newOut['magXarray'][0], newOut['posYarray'][0])
