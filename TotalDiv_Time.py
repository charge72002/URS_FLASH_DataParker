#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  4 17:04:50 2021

@author: Sherry Wong

Code heavily borrowed from EnergyConservation.py and Curl_Div.py
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
sys.path.insert(0, pwd)
# import copy

import athena_read as ath
import h5py
import hdf5_parser

import csv
import time
import os
import beepy
    
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

#%%
startTime = time.time()
timeStamps = [] #list of time stamps
divAvg = [] #list of average divergence
divTotal = [] #list of total divergence
curlAvg = [] #list of average curl
curlTotal = [] #list of total curl

for fileName in sorted(os.listdir(filedirectory)):
    if(fileName.startswith("parkerCRs")):
        #Start - Code from Curl_Div_Gifs.py
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        newOut = hdf5_parser.setup(filedirectory +"/"+ fileName, format = "cartesian")
        x = newOut['posXarray']
        y = newOut['posYarray']
        magY = newOut['magYarray']
        magX = newOut['magXarray']
        
        div = calcDiv(magX, magY)
        curl = calcCurl(magX, magY)
        absDiv = abs(div)
        absCurl = abs(curl)
        
        timeStamps.append(timeStamp)
        divAvg.append(np.mean(absDiv))
        divTotal.append(np.sum(absDiv))
        curlAvg.append(np.mean(absCurl))
        curlTotal.append(np.sum(absCurl))
beepy.beep(6)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Energy conservation done. Time elapsed (sec): " + str(time.time()-startTime))

#%% plot energy total
xticks = np.linspace(0, round(int(max(timeStamps))/10)*10, #round to nearest 10
                     round(int(max(timeStamps))/10)+1)

plt.close()
plt.plot(timeStamps, divAvg)
plt.title("Average Divergence")
plt.xticks(xticks)
plt.ylabel(r"Divergence ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
plt.xlabel("Time stamp")
plt.savefig(pwd + '/Conservation/DivAvg')

plt.close()
plt.plot(timeStamps, curlAvg)
plt.title("Average Curl")
plt.xticks(xticks)
plt.ylabel(r"Curl ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
plt.xlabel("Time stamp")
plt.savefig(pwd + '/Conservation/CurlAvg')

#%% Save as .csv
with open(pwd + '/Conservation/DivCurl.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', "Avg_Divergence", "Total_Divergence", "Avg_Curl", "Total_Curl"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    #iterate through time steps
    for i in range(0, len(timeStamps)):
        writer.writerow({'t':               timeStamps[i], 
                         "Avg_Divergence":  divAvg[i],
                         "Total_Divergence":divTotal[i],
                         "Avg_Curl":        curlAvg[i],
                         "Total_Curl":      curlTotal[i]})
