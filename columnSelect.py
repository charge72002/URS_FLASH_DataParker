# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 23:36:03 2021

@author: wongb
"""

import yt
from yt.units import kpc
import os
from os import path
import numpy as np
import time
import beepy

import moviepy
from moviepy.editor import ImageSequenceClip

os.chdir("/Users/wongb/Documents/Python Scripts")
# filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0020"
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"
# filename = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/parkerCRs_hdf5_plt_cnt_0000"

ds = yt.load(filename)
slc = yt.SlicePlot(ds, 'z', 'density')
slc.save("YT_Test_Plots/HDF5/bad_labels")


for i in sorted(ds.field_list):
    print(i)
for i in sorted(ds.derived_field_list):
    print(i)

###############
# select region
###############
#https://yt-project.org/doc/analyzing/filtering.html#cut-regions
ad = ds.all_data()

bounds = {'xmin': -0.3e+22, 'xmax': 0.22e+22, 'ymin': min(ad['y']).value,'ymax': 0}

dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
dsSelect = dsSelect.cut_region("obj['temp'] > .25*10e3")
slc = yt.SlicePlot(ds, 'z', 'density', data_source=dsSelect, 
                    center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.set_width(-min(ad['y']))

dsSelect = ad.cut_region("obj['temp'] > .35*10e3")
slc = yt.SlicePlot(ds, 'z', 'density', data_source=dsSelect)

# slc.set_ylim('y', min(ad['y']), 0)
slc.annotate_title("Density, selected by temp > $0.35*10^3$")
slc.save("YT_Test_Plots/HDF5/cutRegionE")

#make a lineplot of the region
#https://yt-project.org/doc/visualizing/plots.html#d-line-sampling
slc = yt.LinePlot(ds, 'density', [-0.25e22, 0, 0], [-0.05e22, bounds['ymin'], 0], 512)
slc.save("YT_Test_Plots/HDF5/density_line")

#plot over time
directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
# directory = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/"
# directory = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
bounds = {'xmin': -0.3e+22, 'xmax': 0.22+22, 'ymin': min(ad['y']).value,'ymax': 0}

startTime = time.time()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        ds = yt.load(directory+fileName)
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
        dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
        slc = yt.LinePlot(ds, 'density', [-0.25e22, 0, 0], [-0.05e22, bounds['ymin'], 0], 512)
        if (not path.exists(saveDirectory+"/density_line/")):
            os.mkdir(saveDirectory+"/density_line/")
        slc.save(saveDirectory+"/density_line/" + timeStamp)
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))

##########
# Make Gif
##########
#https://www.tutorialexample.com/python-create-gif-with-images-using-moviepy-a-complete-guide-python-tutorial/

images = []
for fileName in os.listdir(saveDirectory + "/density_line/"): #set in initial parameters
    images.append(saveDirectory + "/density_line/" + fileName)
print(images)
clip = ImageSequenceClip(images, fps=5)
clip.write_gif(saveDirectory + "/density_line/" + '.gif') #saves in outside folder
clip.close()
#os.rmdir(saveDirectory + '/' + field) #delete images to save space