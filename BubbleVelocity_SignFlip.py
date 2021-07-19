# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 16:56:53 2021

@author: wongb
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 19:21:31 2021

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
# import h5py
import numpy as np
# import scipy
from scipy import signal

# import matplotlib
import matplotlib.pyplot as plt #Same installation as above
from matplotlib import ticker, cm #for log scale
import matplotlib.colors as colors #for color mapping

import os
from os import path
import csv

import yt
from yt.data_objects.level_sets.api import Clump, find_clumps
yt.toggle_interactivity()

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import hdf5_parser

#%%
out = hdf5_parser.setup(filename)

ds = yt.load(filename)
ad = ds.all_data()
x = est.closestNum(out['posXarray'], 2.32019*pow(10, 22)) #B
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
ylim = -1.79040182984184e21

#%%
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>-1.0*(10**22)))
ySlice = out['posYarray'][TFselect]
magSlice = out['magXarray'][TFselect]
zeroIndex = est.findSignFlips(magSlice)
magDerivative = np.diff(magSlice, axis=0)
magDerivative = np.insert(magDerivative, 0, 0) #to match shape

#Trim zeroes that are not bubbles
zeroIndexTrimmed = []
for index in zeroIndex:
    if magDerivative[index] > 10**-9: zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)

fig, ax = plt.subplots(2)
ax[0].plot(ySlice, magSlice)
ax[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
for y in ySlice[zeroIndex]:
    ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #prominence lines
ax[1].plot(ySlice, magSlice)
ax[1].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
for y in ySlice[zeroIndexTrimmed]:
    ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #prominence lines
    
#%%
### iterate
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>-1.0*(10**22)))
ySlice = out['posYarray'][TFselect]
for t in range(65, 85):
    filename = filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)
    out = hdf5_parser.setup(filename)
    
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0) #to match shape
    
    zeroIndexTrimmed = []
    for index in zeroIndex:
        if magDerivative[index] > 10**-9: zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    fig, ax = plt.subplots(2)
    fig.suptitle('Zero intercepts, t='+ str(t) +', x='+str(x))
    ax[0].plot(ySlice, magSlice)
    ax[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    for y in ySlice[zeroIndex]:
        ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #prominence lines
    ax[1].plot(ySlice, magSlice)
    ax[1].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    for y in ySlice[zeroIndexTrimmed]:
        ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #prominence lines
    
    fig.savefig(pwd + "/bubble_velocity/col12_t" + str(t))
    plt.close()