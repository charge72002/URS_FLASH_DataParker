# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 13:29:41 2021

@author: Sherry Wong
"""

import yt
from yt.units import kpc
import os
from os import path
import numpy as np
import time
import beepy

import matplotlib
import matplotlib.pyplot as plt

import moviepy
from moviepy.editor import ImageSequenceClip

os.chdir("/Users/wongb/Documents/Python Scripts")
# filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0020"
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"
# filename = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/parkerCRs_hdf5_plt_cnt_0000"

ds = yt.load(filename)
ad = ds.all_data()

for i in sorted(ds.field_list):
    print(i)
for i in sorted(ds.derived_field_list):
    print(i)
    
flux = ad['density']*ad['vely']
plt.clf()
# lev = np.logspace(np.log10(flux.min()), np.log10(flux.max()), num=100)
plt.tricontourf(ad['x'], ad['y'], flux, norm=matplotlib.colors.SymLogNorm(linthresh=0.01))
plt.title("0076: flux using SymLogNorm")
plt.colorbar()
plt.show()
plt.savefig("/Users/wongb/Documents/Python Scripts/YT_Test_Plots/HDF5/0076fluxC")

conversion = 3.086e21
slc = yt.LinePlot(ds, 'vely', [6.25*kpc, -1.75*kpc, 0], [0,  -1.75*kpc, max(ad['x'])*(conversion)], 512)
slc.save("/Users/wongb/Documents/Python Scripts/YT_Test_Plots/HDF5/0076xvely")

# slc = yt.SlicePlot(ds, 'z', flux)
# slc.save("YT_Test_Plots/HDF5/flux")

# for index in range(ad['vely'].size):
#     flux.append(ad['density']*ad['vely'])