#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 08:28:50 2021

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
import matplotlib.colors as colors
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

def EZlabels(title, timestamp, savefile = ""):
    if savefile == "": savefile = title
    plt.colorbar(label = title + r" ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
    plt.title("Mag "+title+" (t="+timestamp+")")
    plt.xlabel("x position (cm)")
    plt.ylabel("y position (cm)")
    plt.savefig(pwd + "/Plots/"+savefile+timestamp+".png")
    
#div
def calcDiv(X_field, Y_field):
    PxFx = np.diff(Y_field, n=1, axis = 0) #shape (511, 512)
    PyFy = np.diff(X_field, n=1, axis = 1) #shape (512, 511)
    return(PxFx[:,1:] + PyFy[1:]) #reselect to rehape

#curl
def calcCurl(X_field, Y_field):
    PxFy = np.diff(Y_field, n=1, axis = 1) #shape (512, 511)
    PyFx = np.diff(X_field, n=1, axis = 0) #shape (511, 512)
    return(PxFy[1:] - PyFx[:,1:]) #reselect to rehape

# Plotting curl and div~~~~~~~~~~~~~~~~~~~~~~
#SETUP DATA

#DISCRETE CURL FORMULA:
#Partial x(Fy) - Partial y(Fx)

timestamp = "0000"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
rawhdf = h5py.File(filename, 'r')
newOut = hdf5_parser.setup(filename, format = "cartesian")

x = newOut['posXarray']
y = newOut['posYarray']
magY = newOut['magYarray']
magX = newOut['magXarray']
magInit = np.sqrt(pow(magX, 2), pow(magY, 2))[1:, 1:]

divInit = calcDiv(magX, magY)
curlInit = curl = calcCurl(magX, magY)


# ~~~~~~~~~~~~~~~~~~~~~~
timestamp = "0080"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
rawhdf = h5py.File(filename, 'r')
newOut = hdf5_parser.setup(filename, format = "cartesian")

x = newOut['posXarray']
y = newOut['posYarray']
magY = newOut['magYarray']
magX = newOut['magXarray']

div = calcDiv(magX, magY)
curl = calcCurl(magX, magY)

#%% Plot (w/ linear scale)
plt.clf()
# lev = np.logspace(np.log10(div.min()), np.log10(div.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.MaxNLocator(100))
EZlabels("Div", timestamp)

plt.clf()
# lev = np.logspace(np.log10(curl.min()), np.log10(curl.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, locator=ticker.MaxNLocator(100))
EZlabels("Curl", timestamp)

#%% LOG SCALE
plt.clf()
div = abs(div)
lowerexp = -11
lowerbound = max(div.min(), pow(10, lowerexp)) #set minimum lower bound
upperbound = max(div.max(), pow(10, lowerexp+1)) #set minimum lower bound
lev = np.logspace(np.log10(lowerbound), np.log10(upperbound), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, levels = lev)
EZlabels("Div", timestamp, savefile="LOGDiv(Exp"+str(lowerexp)+")") 

plt.clf()
curl = abs(curl)
lowerbound = max(curl.min(), pow(10, -9)) #set minimum lower bound
lev = np.logspace(np.log10(lowerbound), np.log10(curl.max()), num=1000)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, levels = lev)
# plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.LogLocator())
EZlabels("Div", timestamp, savefile="LOGCurl(Exp"+str(lowerexp)+")") 

#%% NORMALIZE BY MATH

plt.clf()
z=abs((div-divInit)/magInit)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z), levels=100)
EZlabels("Div", timestamp, savefile="DivSubtractDivideInitLOG")

plt.clf()
z=abs((curl-curlInit)/magInit)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z), levels=100)
EZlabels("Curl", timestamp, savefile="CurlSubtractDivideInitLOG")

#%% experiment with colorbar limits
plt.clf()
fig = plt.figure(1, figsize = (5, 5))
ax = fig.subplots(nrows = 1, ncols = 1)
z=abs((curl-curlInit)/magInit)
#im = ax.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z), levels=100)
# X = x[1:, 1:]
# Y = y[1:, 1:]
X = x[:1, 1:].flatten()
Y = y[1:, :1].flatten()
Z = np.log10(z)
# im = ax.contourf(X,Y, Z, levels=100)
im = ax.contourf(X,Y, Z, levels=100, vmin = -5, vmax = 5)
#ax.imshow(np.log10(z))
# cbar_ax = fig.add_axes([0.91, 0.25, 0.03, 0.5])
# fig.colorbar(im, cax = cbar_ax)
fig.colorbar(im, label = r"log(Curl) ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
# plt.clim(-5, 3)
ax.set_title("Mag Diff(t="+timestamp+")")
ax.set_xlabel("x position (cm)")
ax.set_ylabel("y position (cm)")

plt.show()
#plt.savefig(pwd + "/Plots/Curl_NewColorbar"+timestamp+".png")
#%% PLOT SIDE BY SIDE

#plt.figure(n)

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
    np.reshape(array, (xlen, ylen)) #this throws an error

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

#%% 2D CURL DETECTION
#%%
# plot bubble positions on yt slice
# iteratively test different thresholds
# for threshold in np.linspace(1.0, 7.5, 14):
for threshold in [0, 1, 2, 3]:
#for threshold in np.linspace(100, 500, 1000):
    TFselect = curl > threshold #using this threshold
    print(str(sum(TFselect)) + " matches")
    
    xbubblepos = x[1:,1:][TFselect]
    ybubblepos = y[1:,1:][TFselect]
    
    plt.clf()
    plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), abs((curl-curlInit)/magInit), locator=ticker.MaxNLocator(100))
    EZlabels("Curl", timestamp)
    
    length = 1e21
    #plot '+' on percieved bubble positions
    for i in range(0, len(xbubblepos)):
        x = xbubblepos[i]
        y = ybubblepos[i]
        plt.annotate_line((x, y+length, 0), (x, y-length, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
        plt.annotate_line((x+length, y, 0), (x-length, y, 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    plt.savefig(pwd + "/Plots/2DCurlDetection"+str(threshold)+str(timestamp)+".png")

