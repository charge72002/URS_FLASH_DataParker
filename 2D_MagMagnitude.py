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
yt.toggle_interactivity()

import numpy as np

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
# iteratively test different thresholds
# for threshold in np.linspace(1.0, 7.5, 14):
for threshold in np.linspace(5.0, 10.0, 11):
    TFselect = ad[('gas', 'magnetic_field_magnitude')] < threshold*(10**-9) #using this threshold
    print(str(sum(TFselect)) + " matches")
    
    xbubblepos = ad['x'][TFselect]
    ybubblepos = ad['y'][TFselect]
    
    slc = yt.SlicePlot(ds, 'z', 'magnetic_field_magnitude')
    # slc = yt.SlicePlot(ds, 'z', 'density')
    slc.annotate_streamlines('magnetic_field_x','magnetic_field_y', density=6, plot_args={'linewidth':0.25,'color':'black'})
    length = 1E21*cm
    #plot '+' on percieved bubble positions
    for i in range(0, len(xbubblepos)):
        x = xbubblepos[i]
        y = ybubblepos[i]
        slc.annotate_line((x, y+length, 0), (x, y-length, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
        slc.annotate_line((x+length, y, 0), (x-length, y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    slc.save(pwd + "/bubble_velocity_2D/" + str(threshold))

#%%
## visualize threshold with yt.cut_region
for threshold in np.linspace(1.0, 10.0, 19):
    dsSelect = ad.cut_region("obj['magnetic_field_magnitude'] < " + str(threshold) + "e-9")
    slc = yt.SlicePlot(ds, 'z', 'magnetic_field_x', data_source=dsSelect)
    # slc.annotate_streamlines('magnetic_field_x','magnetic_field_y', density=6, plot_args={'linewidth':0.25,'color':'black'})
    slc.save(pwd + "/bubble_velocity_2D/dsSelect" + str(threshold))

#%%
## what does the magnetic field strength even look like
slc = yt.SlicePlot(ds, 'z', 'magnetic_field_magnitude')
slc.show_colorbar()
slc.save(pwd + "/bubble_velocity_2D/mag")

slc = yt.SlicePlot(ds, 'z', 'magnetic_field_x')
slc.show_colorbar()
slc.save(pwd + "/bubble_velocity_2D/magX00")

slc = yt.SlicePlot(ds, 'z', 'magnetic_field_y')
slc.show_colorbar()
slc.save(pwd + "/bubble_velocity_2D/magY00")

    
#%%
# iterative implementation of threshold
TFselect = ad[('gas', 'magnetic_field_magnitude')] < 4E-9 #using this threshold
print(str(sum(TFselect)) + " matches")
    
slc = yt.SlicePlot(ds, 'z', 'magnetic_field_magnitude')
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y', density=6, plot_args={'linewidth':0.25,'color':'r'})
length = 1E21*cm
num_matches=0
for i in range(0, len(TFselect)):
    if TFselect[i]:
    #select based on mag slope
        avgSlope = 0.5*(ad['magx'][i-1]-ad['magx'][i+1]) #broken. selects Z order neighbors which means nothing
        if avgSlope > 5*(10**-11):
            num_matches+=1
            slc.annotate_line((x, y+length, 0), (x, y-length, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
            slc.annotate_line((x+length, y, 0), (x-length, y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
print(str(sum(TFselect)) + " matches")
print("Reduced to " + str(num_matches) + " matches\n")
slc.save(pwd + "/bubble_velocity_2D/MagMagnitudeTrim")
   