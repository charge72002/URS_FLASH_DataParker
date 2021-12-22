#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 19:02:26 2021

@author: Sherry Wong
"""
#%%
## Windows Filenames
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/wongb/Documents/Python_Scripts"
#%%
## Mac filenames
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0001"
filedirectory = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/bwong/URS_FLASH_DataParker"
#%%

import yt
yt.toggle_interactivity()
import os
import os.path
from os import path
import shutil
import time
import matplotlib
import matplotlib.pyplot as plt
import moviepy
from moviepy.editor import ImageSequenceClip
import beepy #sound for when the code is done running
import numpy as np
from yt.units import kpc
import sys

import sys
sys.path.insert(0, pwd)
import athena_read as ath
import h5py
import hdf5_parser

#%%
timestamp = "0076"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
rawhdf = h5py.File(filename, 'r')
newOut = hdf5_parser.setup(filename, format = "cartesian")

x = newOut['posXarray']
y = newOut['posYarray']
magY = newOut['magYarray']
magX = newOut['magXarray']

#%%
ds = yt.load(filename)
slc = yt.SlicePlot(ds, 'z', 'density')
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=5,plot_args={'linewidth':0.5,'color':'r'}) 
slc.set_zlim('density', 1e-33, 1e-24)
slc.annotate_title(timestamp +" density")

#%%
#follow this link for tutorial:
#https://yt-project.org/doc/cookbook/complex_plots.html
from yt.visualization.api import get_multi_plot #no such function?
fig, axes, colorbars = yt.get_multi_plot(2, 3, colorbar="horizontal", bw=6)
