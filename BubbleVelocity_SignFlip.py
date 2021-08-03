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

def removeDuplicates(arrayIN):
    listOUT = []
    for i in arrayIN:
        if(listOUT.count(i)==0):
            listOUT.append(i)        
    arrayOUT = np.asarray(listOUT)
    return arrayOUT
    #END OF METHOD

#Finds the number in the array closest to the target number.
def closestNum(arrayIN, targetNum, removeDuplicates=True):
    if(removeDuplicates): smallArray = removeDuplicates(arrayIN)
    else: smallArray = arrayIN
    
    currentDev = float("inf")
    result = float("inf")
    #TODO: if you wanna have fun make this a binary search thing
    for i in smallArray:
        newDev = abs(targetNum-i)
        if(newDev < currentDev):
            currentDev =  newDev
            result = i
    if(result == float("inf")):
        raise Exception("No match found")
        
    # #quasi-binary search [[BROKEN]]
    # index = (int) (len(smallArray)/2)
    # while(index > 1):
    #     print(str(index) + ', ' + str(smallArray[index]))
    #     if(targetNum < smallArray[index]):  index = (int) (index * 0.5)
    #     else: index = index = index + (int) ((len(smallArray)-index) * 0.5)
    # print("FINAL: " + str(index) + ', ' + str(smallArray[index]))
    # return result
    #END OF METHOD
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
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<-ylim), (out['posYarray']>ymin))
ySlice = out['posYarray'][TFselect]

#%%
### bubble finding testing
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
### plot promininces?
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
### iterate bubble finding
untrimmed = []
trimmed = []
# trimmedIndices = []

TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>ymin))
ySlice = out['posYarray'][TFselect]
for t in range(65, 85):
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
    
    # # Make a line plot
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
import estimate_bubbleVelocity as est
estVelocities = []
for i in range(0, len(trimmed)):
    currentTimeStep = []
    ypos = ySlice[trimmed[i][0]]
    for j in range(0, len(ypos)):
        y = ypos[j]
        if(abs(y)>3E21): #outside 1kpc of disk
            currentTimeStep.append( est.calcVelocity(y) )
        else: #within 1kpc of disk
            #check if data available for previous timestep
            print("Inside the disk")
            if(i==0): raise Exception("No data for previous v0 and y0!")
            if(j==0): raise Exception("No data for previous y0!")
            currentTimeStep.append( est.calcVelocity(y), y0=ypos[j-1], v0=estVelocities[i-1][1] )
    estVelocities.append(currentTimeStep)
#%%
### save important data as .csv
with open(pwd + '/bubble_velocity/signFlip2.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', 'Ypos', 'magslope', 'expected Yvel']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(0, len(trimmed)):
        writer.writerow({'t': i+65, 
                         'Ypos': ySlice[trimmed[i][0]], 
                         'magslope': trimmed[i][1],
                         'expected Yvel': estVelocities[i]})
        
#%%
### Try a detection little left or right
out = hdf5_parser.setup(filename)
# lev = np.logspace(np.log10(out['densityArray'].min()), np.log10(out['densityArray'].max()), num=1000)
# plt.tricontourf(out['posXarray'], out['posYarray'], out['densityArray'], locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

x = est.closestNum(out['posXarray'], 2.32838E22 - (5*dx)) #2.320198554 is good for col 12
# x = est.closestNum(out['posXarray'], -4.211*(10**21) - (5*dx)) #col 2
# plt.clf()

for t in range(1, 10):
    print("("+ str(t) + ", " + str(x) +")")
    TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>ymin))
    ySlice = out['posYarray'][TFselect]    
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0)
    
    zeroIndexTrimmed = []
    for index in zeroIndex:
        avgSlope = 0.5*(out['magXarray'][index-1]-out['magXarray'][index+1])
        if avgSlope > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    # fig, axs = plt.subplots(2)
    # fig.suptitle('MagX, MagX Derivative, x='+str(x))
    # axs[0].plot(ySlice, magSlice)
    # axs[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1) #zero line
    # axs[1].plot(ySlice, magDerivative)
    # for y in ySlice[zeroIndexTrimmed]:
    #     axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    #     axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
    #     axs[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    #     axs[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
        
    plt.close()
    plt.title('x='+str(x))
    plt.plot(ySlice, magSlice)
    plt.plot([max(ySlice), min(ySlice)], [0, 0], lw=1) #zero line
    for y in ySlice[zeroIndexTrimmed]:
        plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
    plt.show()
    fig.savefig(pwd + "/bubble_velocity/Slice" + str(t))
    
    
    #prep for next loop
    x = est.closestNum(out['posXarray'], x + dx)
    
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

#%%



#####################
# 2D BUBBLE DETECTION
#####################


#Columns are labeled 1-12 left to right up/down/up/down etc
colXbounds = [(0, 0),                       #empty b/c I started at 1
              (-2.41E+22,   -1.97E+22),     #1
              (-2.21E+22,   -1.5E+22),      #2
              (-1.97E+22,   -1.01E+22),     #3
              (-1.11E22,    -6.28E+21),     #4
              (-8.32E+21,   -1.65E+20),     #5
              (-2.69E21,    2E+21),         #6
              (-1.65E+20,   7.16E21),       #7
              (2.62E21,     1.04E+22),      #8
              (8.05E+21,    1.64E+22),      #9
              (7.715E+21,   1.8516E+22),    #10
              (1.29E+22,    2.2E+22),       #11
              (2.1E+22,     2.5E+22)]       #12
#column widths:
# colWidth = [0,
#  4.4e+21,
#  7.01e+21,
#  9.6e+21,
#  4.82e+21,
#  8.155e+21,
#  4.69e+21,
#  7.325e+21,
#  7.78e+21,
#  8.35e+21,
#  1.08e+22,
#  9.1e+21,
#  4e+21]
#%%
### YT annotate lines to visualize where these lineplots are
numlines = 15 #number of lines to the left of x
ylim = -1.79040182984184e21
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])  
slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                   center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
shortX = removeDuplicates(out['posXarray']) #do this ONCE now to save time
x = closestNum(shortX, 2.38*pow(10, 22) - (numlines*dx), removeDuplicates=False) #2.320198554 is good for col 12
# slc.annotate_grids()
# slc.annotate_cell_edges(line_width = 0.00001, alpha = 0.5)

textx = x + (10**21)
for t in range(1, numlines):
    slc.annotate_line((x, -1.8*(10**21), 0), (x, -1.23*(10**22), 0), coord_system="data",  plot_args={"color": "blue"})
    slc.annotate_marker((x, (-2-0.5*t)*(10**21), 0), coord_system="data", plot_args={"color": "red"})
    slc.annotate_text((textx, (-2-0.5*t)*(10**21), 0), str(x), coord_system="data", text_args={"color": "red"})
    x = closestNum(shortX, x + dx, removeDuplicates=False)
slc.save(pwd + "/bubble_velocity_2D/line_annotations")
slc.display()

#%%
#rebuild the lineplot code
out = hdf5_parser.setup(filename)
x = est.closestNum(out['posXarray'], 2.32838*pow(10, 22) - (5*dx)) #2.320198554 is good for col 12
plt.clf()

# TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim))
# ySlice = out['posYarray'][TFselect]    
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>ymin))
ySlice = out['posYarray'][TFselect] #this is in order
magSlice = out['magXarray'][TFselect]
# zeroIndex = est.findSignFlips(magSlice)
# magDerivative = np.diff(magSlice, axis=0)
# magDerivative = np.insert(magDerivative, 0, 0)

# zeroIndexTrimmed = []
# for index in zeroIndex:
#     avgSlope = 0.5*(out['magXarray'][index-1]-out['magXarray'][index+1])
#     if avgSlope > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
# zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)

ySlice.sort()
plt.plot(ySlice, magSlice)
plt.plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
# plt.plot([max(ySlice), min(ySlice)], [0, 0], lw=1) #zero line
# for y in ySlice[zeroIndexTrimmed]:
#     plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
#     plt.plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
# plt.show()    
plt.savefig(pwd + "/bubble_velocity/SliceA")
#%%
### Try a detection little left or right
# lev = np.logspace(np.log10(out['densityArray'].min()), np.log10(out['densityArray'].max()), num=1000)
# plt.tricontourf(out['posXarray'], out['posYarray'], out['densityArray'], locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

x = est.closestNum(out['posXarray'], 2.32838*pow(10, 22) - (5*dx)) #2.320198554 is good for col 12
# x = est.closestNum(out['posXarray'], -4.211*(10**21) - (5*dx)) #col 2
plt.clf()
for t in range(1, 10):
    print("("+ str(t) + ", " + str(x) +")")
    TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>ymax))
    ySlice = out['posYarray'][TFselect]    
    magSlice = out['magXarray'][TFselect]
    zeroIndex = est.findSignFlips(magSlice)
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0)
    
    zeroIndexTrimmed = []
    for index in zeroIndex:
        avgSlope = 0.5*(out['magXarray'][index-1]-out['magXarray'][index+1])
        if avgSlope > 5*(10**-11): zeroIndexTrimmed = np.append(zeroIndexTrimmed, index)
    zeroIndexTrimmed = np.asarray(zeroIndexTrimmed, dtype=int)
    
    fig, axs = plt.subplots(2)
    fig.suptitle('MagX, MagX Derivative, x='+str(x))
    axs[0].plot(ySlice, magSlice)
    axs[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
    axs[1].plot(ySlice, magDerivative)
    for y in ySlice[zeroIndexTrimmed]:
        axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
        axs[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
        axs[1].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1)
        
    fig.savefig(pwd + "/bubble_velocity/Slice" + str(t))
    plt.close()
    
    #prep for next loop
    x = est.closestNum(out['posXarray'], x + dx)