# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 16:56:53 2021

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
dx = out['posXarray'][1]-out['posXarray'][0]
# x = est.closestNum(out['posXarray'], 2.32019*pow(10, 22)) #B
x = est.closestNum(out['posXarray'], 2.2871316*pow(10, 22)) #moved closer to middle
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
ylim = -1.79040182984184e21

bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']>-ylim), (out['posYarray']<ymax))
ySlice = out['posYarray'][TFselect]

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
    if magDerivative[index] > 10**-10: zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
print("untrimmed: " + str(zeroIndex))
print("derivatives: " + str(magDerivative[zeroIndex]))
print("trimmed: " + str(zeroIndexTrimmed))
print("derivatives: " + str(magDerivative[zeroIndexTrimmed]))

#%%
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
untrimmed = []
trimmed = []
# trimmedIndices = []

TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>ymin))
ySlice = out['posYarray'][TFselect]
for t in range(79, 81):
    filename = filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)
    out = hdf5_parser.setup(filename)
    
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    # https://numpy.org/doc/stable/reference/generated/numpy.diff.html
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0) #to match shape
        
    zeroIndexTrimmed = []
    for index in zeroIndex:
        if magDerivative[index] > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    print("~~~ t=" + str(t) + " ~~~")
    print("untrimmed: " + str(zeroIndex))
    print("derivatives: " + str(magDerivative[zeroIndex]))
    print("trimmed: " + str(zeroIndexTrimmed))
    print("derivatives: " + str(magDerivative[zeroIndexTrimmed]))
    untrimmed.append([zeroIndex, magDerivative[zeroIndex]])
    trimmed.append([zeroIndexTrimmed, magDerivative[zeroIndexTrimmed]])
    
    ## Make a line plot
    # fig, ax = plt.subplots(2)
    # fig.suptitle('Zero intercepts, t='+ str(t) +', x='+str(x))
    # ax[0].plot(ySlice, magSlice)
    # ax[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    # for y in ySlice[zeroIndex]:
    #     ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    #     ax[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) 
    # ax[1].plot(ySlice, magSlice)
    # ax[1].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    # for y in ySlice[zeroIndexTrimmed]:
    #     ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    #     ax[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) 
    
    # fig.savefig(pwd + "/bubble_velocity/col12_t" + str(t))
    # plt.close()
    
    # #plot slice in yt
    # ds = yt.load(filename)
    # ad = ds.all_data()
    # dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
    # dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
    # z=out['tempArray']
    # lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
    # slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
    #                         center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
    # slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=6,plot_args={'linewidth':0.5,'color':'r'}) 
    # slc.annotate_title(str(t) +" "+ 'temp')
    # slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
    # slc.set_zlim('temp', 1e2, 1e6)
    # slc.annotate_title("Temp (\N{DEGREE SIGN}K) + bubble positions, t=" + str(t) + "x=" + str(x))

    # #plot bubble position
    # zeroIndexTrimmed = []
    # for index in zeroIndex:
    #     if magDerivative[index] > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    # zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    # for y in ySlice[zeroIndexTrimmed]:
    #     slc.annotate_line((x, bounds['ymin'], 0), (x, bounds['ymax'], 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    #     slc.annotate_line((bounds['xmin'], y, 0), (bounds['xmax'], y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    # slc.save(pwd + "/bubble_velocity/t=" + str(t))
   
print(trimmed)
print(str(trimmed))

#%%
### estimate velocities
estVelocities = []
for i in range(0, len(trimmed)):
    currentTimeStep = []
    for ypos in ySlice[trimmed[i][0]]:
        currentTimeStep.append( est.calcVelocity(ypos) )
    estVelocities.append(currentTimeStep)
#%%
### save important data as .csv
with open(pwd + '/bubble_velocity/signFlip.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', 'Ypos', 'magslope', 'expected Yvel']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(0, len(trimmed)):
        writer.writerow({'t': i+65, 
                         'Ypos': ySlice[trimmed[i][0]], 
                         'magslope': trimmed[i][1],
                         'expected Yvel': estVelocities[i]})
        
#%%
### Try a little left or right
lev = np.logspace(np.log10(out['densityArray'].min()), np.log10(out['densityArray'].max()), num=1000)
plt.tricontourf(out['posXarray'], out['posYarray'], out['densityArray'], locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

# x = est.closestNum(out['posXarray'], 2.32838*pow(10, 22)) - (5*dx) #2.320198554 is good for col 12
x = est.closestNum(out['posXarray'], -4.211*(10**21) - (5*dx)) #col 2
plt.clf()
for t in range(1, 10):
    TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']>-ylim), (out['posYarray']<ymax))
    ySlice = out['posYarray'][TFselect]
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    # magDerivative = np.diff(magSlice, axis=0)
    # magDerivative = np.insert(magDerivative, 0, 0) #to match shape
    
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0)
    
    zeroIndexTrimmed = []
    for index in zeroIndex:
        if magDerivative[index] > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    
    plt.plot(ySlice, magSlice)
    plt.plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    for y in ySlice[zeroIndexTrimmed]:
        plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
    plt.savefig(pwd + "/bubble_velocity/Slice" + str(t))
    plt.close()
    x = est.closestNum(out['posXarray'], x + dx)
    print("("+ str(t) + ", " + str(x) +")")
    
#%%
### YT plot position onto slice
# select col 12
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
# ad = ds.all_data()
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
z=out['tempArray']
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                       center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=6,plot_args={'linewidth':0.5,'color':'r'}) 
# slc.annotate_title(timeStamp +" "+ field)
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
slc.set_zlim('temp', 1e2, 1e6)
# ann_slc = plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev)
# plt.title("Temp (\N{DEGREE SIGN}K)")
for y in ySlice[zeroIndexTrimmed]:
    slc.annotate_line((x, bounds['ymin'], 0), (x, bounds['ymax'], 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    slc.annotate_line((bounds['xmin'], y, 0), (bounds['xmax'], y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
slc.save(pwd + "/bubble_velocity/xPositions")

fig = slc.plots['temp'].figure
# figure = fig.figure
# axes = fig.axes
# colorbar_axes = fig.cax

# fig.show()
# slc.display()

#%%
### iterate yt bubble position on slices
ds = yt.load(filename)
ad = ds.all_data()
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']>-ylim), (out['posYarray']<ymax))
ySlice = out['posYarray'][TFselect]

z=out['tempArray']
for t in range(75, 85):
    filename = filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)
    #find bubble position
    out = hdf5_parser.setup(filename)
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    # magDerivative = np.diff(magSlice, axis=0)
    # magDerivative = np.insert(magDerivative, 0, 0) #to match shape
    
    #plot in yt
    ds = yt.load(filename)
    ad = ds.all_data()
    dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
    dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
    lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
    slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                           center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
    slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=6,plot_args={'linewidth':0.5,'color':'r'}) 
    slc.annotate_title(t +" "+ 'temp')
    slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
    slc.set_zlim('temp', 1e2, 1e6)
    slc.annotate_title("Temp (\N{DEGREE SIGN}K) + bubble positions, t=" + str(t) + "x=" + str(x))

    #plot bubble position
    zeroIndexTrimmed = []
    for index in zeroIndex:
        if magDerivative[index] > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    for y in ySlice[zeroIndexTrimmed]:
        slc.annotate_line((x, bounds['ymin'], 0), (x, bounds['ymax'], 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
        slc.annotate_line((bounds['xmin'], y, 0), (bounds['xmax'], y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    slc.save(pwd + "/bubble_velocity/t=" + str(t))

# fig = slc.plots['temp'].figure
# figure = fig.figure
# axes = fig.axes
# colorbar_axes = fig.cax

# fig.show()
# slc.display()