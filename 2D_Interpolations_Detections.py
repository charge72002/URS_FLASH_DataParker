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

#%% Plotting curl and div
#DISCRETE CURL FORMULA:
#Partial x(Fy) - Partial y(Fx)
timestamp = "0001"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
rawhdf = h5py.File(filename, 'r')
newOut = hdf5_parser.setup(filename, format = "cartesian")

x = newOut['posXarray']
y = newOut['posYarray']
magY = newOut['magYarray']
magX = newOut['magXarray']

# #curl
# PxFy = np.diff(magY, n=1, axis = 0) #shape (511, 512)
# PyFx = np.diff(magX, n=1, axis = 1) #shape (512, 511)
# curl = PxFy[:,1:] - PyFx[1:] #reselect to rehape

# #div
# PyFy = np.diff(magY, n=1, axis = 1) #shape (512, 511)
# PxFx = np.diff(magX, n=1, axis = 0) #shape (511, 512)
# div = PxFx[:,1:] + PyFy[1:]

#div
PxFx = np.diff(magY, n=1, axis = 0) #shape (511, 512)
PyFy = np.diff(magX, n=1, axis = 1) #shape (512, 511)
div = PxFx[:,1:] + PyFy[1:] #reselect to rehape

#curl
PxFy = np.diff(magY, n=1, axis = 1) #shape (512, 511)
PyFx = np.diff(magX, n=1, axis = 0) #shape (511, 512)
curl = PxFy[1:] - PyFx[:,1:]

#%%
plt.clf()
# lev = np.logspace(np.log10(div.min()), np.log10(div.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.MaxNLocator(100))
plt.colorbar()
plt.title("Mag Divergence (t="+timestamp+")")
plt.savefig(pwd + "/Plots/Div"+timestamp+".png")

plt.clf()
# lev = np.logspace(np.log10(curl.min()), np.log10(curl.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, locator=ticker.MaxNLocator(100))
plt.colorbar()
plt.title("Mag Curl (t="+timestamp+")")
plt.savefig(pwd + "/Plots/Curl"+timestamp+".png")

#%% LOG SCALE
plt.clf()
div = abs(div)
lowerexp = -11
lowerbound = max(div.min(), pow(10, lowerexp)) #set minimum lower bound
upperbound = max(div.max(), pow(10, lowerexp+1)) #set minimum lower bound
lev = np.logspace(np.log10(lowerbound), np.log10(upperbound), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, levels = lev)
plt.colorbar(label = r"Divergence ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
plt.title("Mag Divergence (t="+timestamp+")")
plt.xlabel("x position (cm)")
plt.ylabel("y position (cm)")
plt.savefig(pwd + "/Plots/LOGDiv"+timestamp+str(lowerexp)+".png")

plt.clf()
curl = abs(curl)
lowerbound = max(curl.min(), pow(10, -9)) #set minimum lower bound
lev = np.logspace(np.log10(lowerbound), np.log10(curl.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, levels = lev)
# plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.LogLocator())
plt.colorbar(label = r"Curl ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
plt.title("Mag Curl (t="+timestamp+")")
plt.xlabel("x position (cm)")
plt.ylabel("y position (cm)")
plt.savefig(pwd + "/Plots/LOGCurl"+timestamp+".png")

#%% Zoom (Why? To re-draw the colorbar while plotting, NOT just matplotlib zoom.)
#trim to length
xtrim = x[1:, 1:]
ytrim = y[1:, 1:]
#set box bounds
conversion = 3.086e21 #(pc to cm?)
ylim = -1.79040182984184e21
bounds = {'xmin': 2.5*conversion, 'xmax': 6*conversion, 'ymin': np.min(y),'ymax': ylim}
#apply bounds
TFtable = np.logical_and(
    np.logical_and(bounds['xmin']<xtrim, bounds['xmax']>xtrim),
    np.logical_and(bounds['ymin']<ytrim, bounds['ymax']>ytrim)
    )
#shape is not preserved in the following operation
xZoom = xtrim[TFtable]
yZoom = ytrim[TFtable]
curlZoom = curl[TFtable]
divZoom = div[TFtable]
#reshape into 2D; need to find box dimensions
middleXpos=(bounds['xmin']+bounds['xmax'])/2
xindex = np.where(xtrim==closestNum(xtrim, middleXpos))
yindex = closestNum(ytrim, (bounds['ymin']+bounds['ymax'])/2)
xlen = np.count_nonzero(TFtable[xindex])
ylen = np.count_nonzero(TFtable[:, yindex])
for array in [xZoom, yZoom, curlZoom, divZoom]:
    np.reshape(array, (xlen, ylen))

#plot
plt.clf()
fig = plt.figure()
ax = plt.contourf(xZoom, yZoom, divZoom, locator=ticker.MaxNLocator(100))
# plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, locator=ticker.MaxNLocator(100))
fig.colorbar(ax, boundaries = np.linspace(-2.5e-7, 2e-7, 10))
plt.title("Mag Divergence (t=80)")
# plt.set_zlim(-3e-7, 3e-7)
plt.savefig("Div_t=80.png")

#%%
# plt.clf()
# plt.imshow(curl) #plots based on index array

# z=magInterp.values.flatten()
# lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.tricontourf(x.flatten(), y.flatten(), z, locator=ticker.LogLocator(), levels=lev)
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

