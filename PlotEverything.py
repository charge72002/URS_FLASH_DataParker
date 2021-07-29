# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 16:38:23 2021

@author: SherryWong
"""

import yt
import os
import os.path
from os import path
import shutil
import time
import matplotlib
import matplotlib.pyplot as plt
import moviepy
from moviepy.editor import ImageSequenceClip
import beepy #sound for when the code is done running
import numpy as np
from yt.units import kpc

# directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
directory = "D://URS_LargeData/Parker_forSherry/"
os.path.exists(directory)
# directory = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/"
# directory = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
path.exists(saveDirectory)


#%%
# # how to get a fixed colorbar???
#slice plotting
# field = "(hrat)"
# field = ('gas', 'hrat')
# field = 'density'
field = 'temp'
startTime = time.time()
#https://stackabuse.com/creating-and-deleting-directories-with-python/
if (not path.exists(saveDirectory + "/" + field)):
    os.mkdir(saveDirectory + "/" + field)
    
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
        
        slc.set_zlim('temp', 1e2, 1e6)
        #Finish
        slc.annotate_title(timeStamp +" "+ field)
        slc.save(saveDirectory + "/" + field + "/" + timeStamp)

        #plot extra per loop        
        # slc.set_zlim(('gas', 'hrat'), 0.1e-27, 1e-28)
        slc = yt.SlicePlot(ds, 'z', 'density')
        slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=5,plot_args={'linewidth':0.5,'color':'r'}) 
        slc.set_zlim('density', 1e-33, 1e-24)
        slc.annotate_title(timeStamp +" density")
        slc.save(saveDirectory + "/density/" + timeStamp)
        
        #plot extra per loop        
        slc = yt.SlicePlot(ds, 'z', 'cray')
        slc.set_zlim('cray', 1e10, 1e13)
        slc.annotate_title(timeStamp +" cray")
        slc.save(saveDirectory + "/cray/" + timeStamp)
        
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))

#%%

#####################
#Calculate total mass
#####################
startTime = time.time()
timeStamps = []
mass = []
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        timeStamps.append(int(timeStamp))
        mass.append(sum(ad[('gas', 'cell_mass')]))
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Mass calc done. Time elapsed (sec): " + str(time.time()-startTime))

plt.clf()
# upperTick = round( np.amax(mass) )
# bottomTick = round( np.amin(mass) )
# log = np.logspace(bottomTick, upperTick, 8)
# plt.plot(timeStamps, np.log10(mass))
plt.plot(timeStamps, mass)
# plt.semilogy(timeStamps, mass)
# plt.set_yscale('log', base=10)


plt.xlabel("Timestamp")
plt.ylabel("Total mass (g)")
# plt.ticklabel_format(axis="y", style="sci", scilimits=(38, 42), useMathText=True)
# plt.ticklabel_format(axis="x", style="plain")
# plt.yticks()
# fileName = "dm1.5_c1_16x16_128x128_Rodrigues_Streaming"
fileName = "m2_c1_16x8_64x64"
plt.savefig(saveDirectory + "/totalMassABC_"+fileName+".png")
plt.show()
#%%
###############################
#Calculate total mass SELECTION
###############################
# pwd = "/Users/bwong/URS_FLASH_DataParker"
pwd = "/Users/wongb/Documents/Python_Scripts"
import sys
sys.path.insert(0, pwd)
# Sherry packages
import hdf5_parser

startTime = time.time()
timeStamps = []
mass = []
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()
ymax = max(ad['y'])
ymin = min(ad['y'])
dsSelect = ad.include_inside('y', ymin, ymax)
# dsSelect = ad.cut_region(["ad['y'] > ymin", "ad['y'] > ymax"])
dsSelect = dsSelect.save_as_dataset(fields=[('gas', 'cell_mass')])


for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        out = hdf5_parser.setup(fileName, [""])
        timeStamps.append(int(timeStamp))
        mass.append(sum(ad[('gas', 'cell_mass')]))
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Mass calc done. Time elapsed (sec): " + str(time.time()-startTime))

plt.clf()
plt.plot(timeStamps, mass)

plt.xlabel("Timestamp")
plt.ylabel("Total mass (g)")
fileName = "m2_c1_16x8_64x64"
plt.savefig(saveDirectory + "/totalMassSelect_m2_c1_16x8_64x64.png")
plt.show()


#%%
##########
# Make Gif
##########
#https://www.tutorialexample.com/python-create-gif-with-images-using-moviepy-a-complete-guide-python-tutorial/

images = []
for fileName in os.listdir(saveDirectory + '/' + field): #set in initial parameters
    images.append(saveDirectory + '/' + field + '/' + fileName)
print(images)
clip = ImageSequenceClip(images, fps=15)
clip.write_gif(saveDirectory + '/' + field + '.gif') #saves in outside folder
clip.close()
#DANGER!!! The line below removes the directory.
# shutil.rmtree(saveDirectory + '/' + field) #delete images to save space

#%%
################
# Zoom CORRECTLY
################

ds = yt.load(directory+"parkerCRs_hdf5_plt_cnt_0000")
ad = ds.all_data()
# field = 'density'
field = 'cray'
tempThreshold = ".35*10e3" #Format as string
conversion = 3.086e21
bounds = {'xmin': 2.5*conversion, 'xmax': 6*conversion, 'ymin': float(min(ad['y']).value),'ymax': 1.79040182984184e21}
ylim = -1.79040182984184e21
# bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}

startTime = time.time()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        
        dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
        dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
        #dsSelect = dsSelect.cut_region("obj['temp'] > " + tempThreshold)
        # dsSelect = dsSelect.cut_region("obj['y'] < 0")
        slc = yt.SlicePlot(ds, 'z', field, data_source=dsSelect, 
                   center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
        slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))

        slc.annotate_title(timeStamp +" "+ field)
        #streamlines
        slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
        # plot = slc.plots[field]
        # slc.set_zlim('density', 1e-33, 1e-24)
        # slc.set_zlim('temp', 1e2, 1e6)
        slc.set_zlim('cray', 1e10, 1e13)
        
        if (not path.exists(saveDirectory + "/" + field)):
            os.mkdir(saveDirectory + "/" + field)
        slc.save(saveDirectory + "/" + field + "/" + timeStamp)
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + str(time.time()-startTime))

#%%
#This actually opens a new display window with YT!!!
#Plot one zoomed image
yt.toggle_interactivity()
#%%
fileName = "parkerCRs_hdf5_plt_cnt_0076"
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"
ds = yt.load(filename)
ad = ds.all_data()
conversion = 3.086e21
bounds = {'xmin': 2.5*conversion, 'xmax': 6*conversion, 'ymin': float(min(ad['y']).value)+0.01e21,'ymax': -1.8e21}
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
slc = yt.SlicePlot(ds, 'z', 'cray', data_source=dsSelect, 
                   center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
# slc.annotate_title("Gas Density, t = 760 Myr")
#slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'})
slc.hide_axes(draw_frame=True)
slc.set_zlim('cray', 1e10, 1e13)
slc.display()
#slc.save(saveDirectory + "/density/" + "test")


# slc = yt.SlicePlot(ds, 'z', 'cray')

slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/0076Cray")
#%%

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
