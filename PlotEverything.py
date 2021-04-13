# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 16:38:23 2021

@author: SherryWong
"""

import yt
import os
import os.path
from os import path
import time
import matplotlib
import matplotlib.pyplot as plt
import moviepy
from moviepy.editor import ImageSequenceClip
import beepy #sound for when the code is done running
import numpy as np

directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
# directory = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/"
# directory = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
path.exists(saveDirectory)
startTime = time.time()

# # how to get a fixed colorbar???
#slice plotting
field = "(hrat)"
field = ('gas', 'hrat')
field = 'density'

for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        
        #SlicePlot
        slc = yt.SlicePlot(ds, 'z', field)
        # slc.annotate_velocity(factor=16)
        slc.annotate_title(timeStamp +" "+ field)
        plot = slc.plots[field]
        slc.set_zlim('density', 1e-33, 1e-24)
        # slc.set_zlim('temp', 1e2, 1e6)
        # slc.set_zlim(('gas', 'hrat'), 0.1e-27, 1e-28)
        # slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 

        #Finish
        #https://stackabuse.com/creating-and-deleting-directories-with-python/
        if (not path.exists(saveDirectory + "/" + field)):
            os.mkdir(saveDirectory + "/" + field)
        slc.save(saveDirectory + "/" + field + "/" + timeStamp)
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))

#####################
#Calculate total mass
#####################
time = []
mass = []
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        time.append(int(timeStamp))
        mass.append(sum(ad[('gas', 'cell_mass')]))
beepy.beep(4)
# upperTick = round( np.amax(mass) )
# bottomTick = round( np.amin(mass) )
# log = np.logspace(bottomTick, upperTick, 8)
plt.plot(time, mass)
plt.xlabel("Timestamp")
plt.ylabel("Total mass (g)")
# plt.ticklabel_format(axis="y", style="sci", scilimits=(38, 42), useMathText=True)
# plt.ticklabel_format(axis="x", style="plain")
# plt.yticks()
fileName = "dm1.5_c1_16x16_128x128_Rodrigues_Streaming"
plt.savefig(saveDirectory + "/totalMass_"+fileName+".png")


##########
# Make Gif
##########
#https://www.tutorialexample.com/python-create-gif-with-images-using-moviepy-a-complete-guide-python-tutorial/

images = []
for fileName in os.listdir(saveDirectory + '/' + field): #set in initial parameters
    images.append(saveDirectory + '/' + field + '/' + fileName)
print(images)
clip = ImageSequenceClip(images, fps=5)
clip.write_gif(saveDirectory + '/' + field + '.gif') #saves in outside folder
clip.close()
#os.rmdir(saveDirectory + '/' + field) #delete images to save space

################
# Zoom CORRECTLY
################

field = 'density'

startTime = time.time()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        
        bounds = {'xmin': 2.75e+22, 'xmax': 4.5+22, 'ymin': min(ad['y']).value,'ymax': 0}
        dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
        dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
        # dsSelect = dsSelect.cut_region("obj['y'] < 0")
        slc = yt.SlicePlot(ds, 'z', field, data_source=dsSelect, 
                   center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
        slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))

        slc.annotate_title(timeStamp +" "+ field)
        plot = slc.plots[field]
        slc.set_zlim('density', 1e-33, 1e-24)
        
        if (not path.exists(saveDirectory + "/" + field)):
            os.mkdir(saveDirectory + "/" + field)
        slc.save(saveDirectory + "/" + field + "/" + timeStamp)
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))



# =============================================================================
# timeStamp = "0076"
# ds = yt.load(directory+fileName)
# slc = yt.SlicePlot(ds, 'z', field)
# slc.annotate_title(timeStamp +" "+ field)
# plot = slc.plots[field]
# plot.cb.set_ticks([1e3, 1e4, 1e5])
# plot.cb.set_ticklabels(["$10^3$", "$10^4$", "$10^5$"])
# slc.save("YT_Test_Plots/HDF5/ColorbarTest")
# # slc.save(saveDirectory + "/" + field + "/" + timeStamp)
# 
# mass = []
# # Cell Mass
# for fileName in os.listdir(directory):
#     if(fileName.startswith("parkerCRs")):
#         #Start
#         print(fileName)
#         ds = yt.load(directory+fileName)
#         print(sum(ds[('gas', 'cell_mass')]))
# =============================================================================
