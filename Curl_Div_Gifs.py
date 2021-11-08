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


#%%
startTime = time.time()
#https://stackabuse.com/creating-and-deleting-directories-with-python/
    
for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        newOut = hdf5_parser.setup(filedirectory +"/"+ fileName, format = "cartesian")
        x = newOut['posXarray']
        y = newOut['posYarray']
        magY = newOut['magYarray']
        magX = newOut['magXarray']
        
        PxFy = np.diff(magY, n=1, axis = 0) #shape (511, 512)
        PyFx = np.diff(magX, n=1, axis = 1) #shape (512, 511)
        curl = PxFy[:,1:] - PyFx[1:] #reselect to rehape
        
        #div
        PyFy = np.diff(magY, n=1, axis = 1) #shape (512, 511)
        PxFx = np.diff(magX, n=1, axis = 0) #shape (511, 512)
        div = PxFx[:,1:] + PyFy[1:]
        
        plt.clf()
        fig = plt.figure() 
        ax = plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), div, locator=ticker.MaxNLocator(100))
        plt.colorbar(ax)
        plt.title("Mag Divergence (t="+str(timeStamp)+")")
        # ax.set_zlim(-3e-7, 3e-7)
        plt.savefig(saveDirectory + "/Div/t="+str(timeStamp)+".png")
        
        plt.clf()
        fig = plt.figure() 
        plt.contourf(x[:1, 1:].flatten(), y[1:, :1].flatten(), curl, locator=ticker.MaxNLocator(100))
        plt.colorbar()
        plt.title("Mag Curl (t="+str(timeStamp)+")")
        # plt.set_zlim(-2.5e-7, 2.5e-7) #
        # set_zlim doesn't really work; 
        # set colobar boundaries/contour levels manually?
        plt.savefig(saveDirectory + "/Curl/t="+str(timeStamp)+".png")
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))

#%%
##########
# Make Gif
##########
#https://www.tutorialexample.com/python-create-gif-with-images-using-moviepy-a-complete-guide-python-tutorial/

field = "Curl"
images = []
for fileName in os.listdir(saveDirectory + '/' + field): #set in initial parameters
    if (fileName.endswith(".png")):#avoid file format errors, even hidden 
        images.append(saveDirectory + '/' + field + '/' + fileName)
#USE GLOB TO ORDER FILES
print(images)
clip = ImageSequenceClip(images, fps=15)
clip.write_gif(saveDirectory + '/' + field + '.gif') #saves in outside folder
clip.close()
#DANGER!!! The line below removes the directory.
# shutil.rmtree(saveDirectory + '/' + field) #delete images to save space


