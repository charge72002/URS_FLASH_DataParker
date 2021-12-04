#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 18:53:49 2021

@author: Sherry Wong

This code comes from:
    "2D_Interpolations_Detections.py"
    "PlotEverything.py"
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
saveDirectory = pwd + "/gifs"

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

import os
from os import path
import time
import beepy #sound for when the code is done running
import moviepy
from moviepy.editor import ImageSequenceClip
import shutil
# import glob #order files


def directoryCheck(testDir):
    #make directory if they do not yet already exist
    if (not path.exists(testDir)):
        os.mkdir(testDir)
        
    #ask if usr wants to overwrite existing files
    if(len(os.listdir(testDir)) != 0):
        while(True):
            print("Files found in " + testDir)
            print("Overwrite? (Y/N)", end="")
            response = input()
            print() #linebreak
            if response == 'Y': break
            if response == 'N': 
                print("Ending program. Reason: Files found in " + testDir)
                sys.exit()
            elif response != 'Y': print("Invalid response.")
        print("Emptying "+ testDir)
        #DANGER!!! The line below removes the directory.
        shutil.rmtree(testDir) #remove path
        os.mkdir(testDir) #remake path

#taken from Curl_Div.py        
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
curlInit = calcCurl(magX, magY)

#%%
while(True):
    print("Files found in " + saveDirectory + "/curl")
    print("Overwrite? (Y/N)", end="")
    response = input()
    print() #linebreak
    if response == 'Y': break
    if response == 'N': 
        print("Ending program. Reason: Files found in " + saveDirectory + "/curl")
        sys.exit()
    elif response != 'Y': print("Invalid response.")
print("Emptying "+ saveDirectory + "/curl")
#%%
startTime = time.time()
#https://stackabuse.com/creating-and-deleting-directories-with-python/

directoryCheck(saveDirectory + "/curl") 
directoryCheck(saveDirectory + "/div") 
    while(True):
        print("Files found in " + testDir)
        print("Overwrite? (Y/N)", end="")
        response = input()
        print() #linebreak
        if response == 'Y': break
        if response == 'N': 
            print("Ending program. Reason: Files found in " + testDir)
            sys.exit()
        elif response != 'Y': print("Invalid response.")
    print("Emptying "+ testDir)
    #DANGER!!! The line below removes the directory.
    shutil.rmtree(testDir) #remove path
    os.mkdir(testDir) #remake path
    
for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):            
        plt.figure(1) #limit to one plot window
        
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        newOut = hdf5_parser.setup(filedirectory +"/"+ fileName, format = "cartesian")
        x = newOut['posXarray']
        y = newOut['posYarray']
        magY = newOut['magYarray']
        magX = newOut['magXarray']
        
        div = calcDiv(magX, magY)
        curl = calcCurl(magX, magY)
        
        plt.clf()
        # fig = plt.figure()
        # ax = fig.add_subplot(1, 1, 1)

        ##CURL
        z = abs((curl-curlInit)/magInit)
        #lowerbound = max(z.min(), pow(10, -9)) #set minimum lower bound
        lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=100)
        ax = plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z),levels=100)#, levels = lev)
        # plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.LogLocator())
        plt.title("Mag Curl Diff/Initial Mag Field"+"(t="+str(timeStamp)+")")
        plt.xlabel("x position (cm)")
        plt.ylabel("y position (cm)")
        # ax.set_zlim(-9, 5)
        plt.colorbar(label =  r"log(Curl) ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
        plt.savefig(saveDirectory + "/curl/t="+str(timeStamp)+".png")
        
        ##DIV
        z = abs((div-divInit)/magInit)
        #lowerbound = max(z.min(), pow(10, -9)) #set minimum lower bound
        lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=100)
        ax = plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), np.log10(z),levels=100)#, levels = lev)
        # plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.LogLocator())
        plt.title("Mag Div Diff/Initial Mag Field"+"(t="+str(timeStamp)+")")
        plt.xlabel("x position (cm)")
        plt.ylabel("y position (cm)")
        # ax.set_zlim(-9, 5)
        plt.colorbar(label =  r"log(Div) ($\frac{G}{cm} = \frac{1}{10^{-7}} \frac{\mu G}{kpc}$)")
        plt.savefig(saveDirectory + "/div/t="+str(timeStamp)+".png")
        
        plt.clf() #close plot window
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))
plt.close()

#%%
##########
# Make Gif
##########
#https://www.tutorialexample.com/python-create-gif-with-images-using-moviepy-a-complete-guide-python-tutorial/

field = "div"
images = []
for fileName in os.listdir(saveDirectory + '/' + field): #set in initial parameters
    if (fileName.endswith(".png")):#avoid file format errors, even hidden 
        images.append(saveDirectory + '/' + field + '/' + fileName)
#USE GLOB TO ORDER FILES
images = sorted(images)
print(images)
clip = ImageSequenceClip(images, fps=15)
# clip.write_gif(saveDirectory + '/' + field + '.gif') #saves in outside folder
clip.write_videofile(saveDirectory + '/' + field + '.mp4') #save as .mp4

clip.close()
#DANGER!!! The line below removes the directory.
# shutil.rmtree(saveDirectory + '/' + field) #delete images to save space



