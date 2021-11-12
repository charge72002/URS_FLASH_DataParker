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
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0001"
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
import sympy

#Methods from "2D_xy_Scans.py"
##THIS HAS BEEN ADAPTED FOR 2DARRAYS
def removeDuplicates(arrayIN):
    listOUT = []
    for j in arrayIN:
        for i in j:
            if(listOUT.count(i)==0):
                listOUT.append(i)        
    arrayOUT = np.asarray(listOUT)
    return arrayOUT
    #END OF METHOD
    
#Finds the number in the array closest to the target number.
#Requires removeDuplicates() function
##THIS HAS BEEN ADAPTED FOR 2DARRAYS
def closestNum(arrayIN, targetNum, removeDupe=True):
    arrayIN = list(arrayIN) #remove casting errors
    if(removeDupe): smallArray = removeDuplicates(arrayIN)
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

#%%
rawhdf = h5py.File(filename, 'r')
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

#%% 1D interpolate

#row 0
row = 200
x = newOut['posXarray'][row]
magY = newOut['magYarray'][row]
Xslice = interp.interp1d(x, magY)

# y = newOut['posYarray'][0]
# magX = newOut['magXarray'][0]
# Yslice = interp.interp1d(y, magX)

xnew = np.linspace(x.min(), x.max(), 512)
ynew = Xslice(xnew)   # use interpolation function returned by `interp1d`
plt.clf()
plt.plot(x, magY, 'o', xnew, ynew, '-')
plt.plot([x.min(), x.max()], [0, 0])
# plt.plot(x, magY)
plt.savefig(pwd + "/Plots1Dinterp.png")


#%%
#plot where the lines are
x = newOut['posXarray'].flatten()
y = newOut['posYarray'].flatten()
z = newOut['densityArray'].flatten()
plt.clf()
plt.tricontourf(x, y, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
# plt.plot([x.min(), x.max()], [y[row][0], y[row][0]])
plt.savefig(pwd + "/Plots/Xslice_t=80.png")


#%% 2D interpolate
## find zeros
x = newOut['posXarray']
y = newOut['posYarray']
magY = newOut['magYarray']
magX = newOut['magXarray']
pow(np.asarray([1, 2, 3]), 2)
mag = np.sqrt(pow(magX, 2) + pow(magY, 2))
# magInterp = interp.RectBivariateSpline(x, y, mag)
# magInterp = interp.interp2d(x, y, mag, kind="linear")

#select list of x values, list of y values
magInterp = interp.RegularGridInterpolator((x[0, :], y[:,0]), mag) #this works!
#CURL AND DIVERGENCE
grad = np.gradient(magInterp)

#%% TEST: Plotting interpolated vs. original data (sanity check)
#Compare interpolated and original data plots. They should match.
#This code takes a few moments (<1min) to run
plt.clf()
fig, axs = plt.subplots(2, 1) #print on the same axes
fig.suptitle("Mag magnitude")

axs[0].set_title("Interpolated Data")
interp=magInterp.values.flatten() #must be flattened b/c tricountour prefers 1D data
lev = np.logspace(np.log10(interp.min()), np.log10(interp.max()), num=1000)
axs[0].tricontourf(x.flatten(), y.flatten(), interp, locator=ticker.LogLocator(), levels=lev)

axs[1].set_title("Original Data")
z = mag.flatten()
levz = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
axs[1].tricontourf(x.flatten(), y.flatten(), z, locator=ticker.LogLocator(), levels=levz)

plt.savefig(pwd+"/Plots/InterpCheck.png")
## check neighbors


#%% check increasing/decreasing x, y
x = newOut['posXarray'].flatten()
for index in range(1, len(x)):
    if(x[index-1] >= x[index]):
        print("ERROR")
        print(index-1)
        print(x[index-1])
        print(index)
        print(x[index])
        break;

