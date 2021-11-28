#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 20:12:23 2021

@author: Sherry Wong

This code applies the concepts from "2D_MagMagnitude.py" to the "Curl_Div.py" plots.
The goal is to detect bubbles (magnetic o points) as they evolve in the simulation.
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
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import h5py
import hdf5_parser

#%%
#div
#div
def calcDiv(X_field, Y_field):
    PxFx = np.diff(magY, n=1, axis = 0) #shape (511, 512)
    PyFy = np.diff(magX, n=1, axis = 1) #shape (512, 511)
    return(PxFx[:,1:] + PyFy[1:]) #reselect to rehape

#curl
def calcCurl(X_field, Y_field):
    PxFy = np.diff(magY, n=1, axis = 1) #shape (512, 511)
    PyFx = np.diff(magX, n=1, axis = 0) #shape (511, 512)
    return(PxFy[1:] - PyFx[:,1:]) #reselect to rehape

def EZlabels(title, timestamp, savefile = ""):
    if savefile == "": savefile = title
    plt.colorbar(label = title + r" ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
    plt.title("Mag "+title+" (t="+timestamp+")")
    plt.xlabel("x position (cm)")
    plt.ylabel("y position (cm)")
    plt.savefig(pwd + "/Plots/"+savefile+timestamp+".png")
#%% Calculate inital curl div

#import data
timestamp = "0000"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
newOut = hdf5_parser.setup(filename, format = "cartesian")
magY = newOut['magYarray']
magX = newOut['magXarray']
magInit = np.sqrt(pow(magX, 2), pow(magY, 2))[1:, 1:]

divInit = calcDiv(magX, magY)
curlInit = calcCurl(magX, magY)

#load coordinates; you only need to do this once
x = newOut['posXarray']
y = newOut['posYarray']
#%% Select time stamp and calculate curl div

#import data
timestamp = "0080"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_" + timestamp
newOut = hdf5_parser.setup(filename, format = "cartesian")
magY = newOut['magYarray']
magX = newOut['magXarray']
magInit = np.sqrt(pow(magX, 2), pow(magY, 2))[1:, 1:]

divNow = calcDiv(magX, magY)
curlNow = calcCurl(magX, magY)

#%% Plot it!
# NORMALIZE BY MATH

# plt.clf()
# z = abs((divNow-divInit)/magInit)
# plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z),levels=100)
# EZlabels("Div", timestamp, savefile="DivSubtractDivideInitLOG")

plt.clf()
z = abs((curlNow-curlInit)/magInit)
z = np.log10(z)
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), z,levels=100)
EZlabels("LogCurl", timestamp, savefile="CurlSubtractDivideInitLOG")
# plt.colorbar(label = r" ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")

#%% Bubble detection

#in 10^threshold
threshold = 2.5
TFselect = (z > threshold) #using this threshold
print(True in TFselect)
print(str(np.sum(TFselect)) + " matches")
xbubblepos = x[1:,1:][TFselect]
ybubblepos = y[1:,1:][TFselect]

# Plot identified bubbles

plt.close()
plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), z,levels=100)
length = .1e22
for coordIndex in range(0, len(xbubblepos)):
    X = xbubblepos[coordIndex]
    Y = ybubblepos[coordIndex]
    # plt.plot((X, Y+length, 0), (X, Y-length, 0),  color="red", linewidth= 1)
    # plt.plot((X+length, Y, 0), (X-length, Y, 0),  color="blue", linewidth= 1)
    plt.plot([X+length, X-length], [Y, Y],  color="red", linewidth= 1)
    plt.plot([X, X], [Y+length, Y-length],  color="red", linewidth= 1)
plt.colorbar(label = r"Log Curl ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
plt.title("Mag Curl Detect, thresh="+str(threshold)+" (t="+timestamp+")")
plt.xlabel("x position (cm)")
plt.ylabel("y position (cm)")
plt.savefig(pwd+"/curl_2D/t=80_thresh="+str(threshold)+".png")

#%% Try it in a gif

import os
from os import path
import time
import beepy #sound for when the code is done running
import moviepy
from moviepy.editor import ImageSequenceClip

for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        newOut = hdf5_parser.setup(filedirectory +"/"+ fileName, format = "cartesian")
        x = newOut['posXarray']
        y = newOut['posYarray']
        magY = newOut['magYarray']
        magX = newOut['magXarray']
        
        z = np.log10(abs((curlNow-curlInit)/magInit))