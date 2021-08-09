# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 14:47:52 2021

@author: wongb
"""
#%%
## Windows Filenames
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
pwd = "/Users/wongb/Documents/Python_Scripts"

#%%
## Mac filenames
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/"
pwd = "/Users/bwong/URS_FLASH_DataParker"

#%%
import yt
import matplotlib.pyplot as plt
from yt.units import cm

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import hdf5_parser

#%%

ds = yt.load(filename)
ad = ds.all_data()
out = hdf5_parser.setup(filename)
# ymax=float(max(ad['y']).value)
# ymin=float(min(ad['y']).value)
# xmax=float(max(ad['x']).value)
# xmin=float(min(ad['x']).value)

#%%
plt.clf()
plt.plot(ad[('gas', 'magnetic_field_magnitude')])
# plt.semilogy(ad[('gas', 'magnetic_field_magnitude')])
plt.plot([0, len(ad[('gas', 'magnetic_field_magnitude')])], [0, 0], lw=1) #lines [x1, x2] [y1, y2]
#select only bubbles
# TFselect = ad[('gas', 'magnetic_field_magnitude')] < 0.25E-9 #using this threshold

#%%
# plot bubble positions on yt slice
TFselect = ad[('gas', 'magnetic_field_magnitude')] < 0.5E-9 #using this threshold
print(str(sum(TFselect)) + " matches")

xbubblepos = ad['x'][TFselect]
ybubblepos = ad['y'][TFselect]

slc = yt.SlicePlot(ds, 'z', 'magnetic_field_magnitude')
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y', density=6, plot_args={'linewidth':0.25,'color':'r'})
length = 1E21*cm
for i in range(0, len(xbubblepos)):
    x = xbubblepos[i]
    y = ybubblepos[i]
    slc.annotate_line((x, y+length, 0), (x, y-length, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    slc.annotate_line((x+length, y, 0), (x-length, y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
slc.save(pwd + "/bubble_velocity_2D/MagMagnitude")

#%%
# iterative implementation
TFselect = ad[('gas', 'magnetic_field_magnitude')] < 0.25E-9 #using this threshold
print(str(sum(TFselect)) + " matches")

slc = yt.SlicePlot(ds, 'z', 'magnetic_field_magnitude')
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y', density=6, plot_args={'linewidth':0.25,'color':'r'})
length = 1E21*cm
for i in range(0, len(TFselect)):
    if TFselect[i]:
    #select based on mag slope
        avgSlope = 0.5*(ad['magx'][i-1]-ad['magx'][i+1])
        if avgSlope > 5*(10**-11):
            slc.annotate_line((x, y+length, 0), (x, y-length, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
            slc.annotate_line((x+length, y, 0), (x-length, y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
slc.save(pwd + "/bubble_velocity_2D/MagMagnitude")
   