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

directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
path.exists(saveDirectory)
startTime = time.time()

# # how to get a fixed colorbar???
#slice plotting
field = "density"
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
        # plot.cb.set_ticks([1e3, 1e4, 1e5])
        # plot.cb.set_ticklabels(["$10^3$", "$10^4$", "$10^5$"])
        # slc.set_zlim('temp', 1e2, 1e6)
        slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 

        #Finish
        #https://stackabuse.com/creating-and-deleting-directories-with-python/
        if (not path.exists(saveDirectory + "/" + field)):
            os.mkdir(saveDirectory + "/" + field)
        slc.save(saveDirectory + "/" + field + "/" + timeStamp)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting done. Time elapsed (sec): " + (time.time()-startTime))

timeStamp = "0076"
ds = yt.load(directory+fileName)
slc = yt.SlicePlot(ds, 'z', field)
slc.annotate_title(timeStamp +" "+ field)
plot = slc.plots[field]
plot.cb.set_ticks([1e3, 1e4, 1e5])
plot.cb.set_ticklabels(["$10^3$", "$10^4$", "$10^5$"])
slc.save("YT_Test_Plots/HDF5/ColorbarTest")
# slc.save(saveDirectory + "/" + field + "/" + timeStamp)

mass = []
# Cell Mass
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        ds = yt.load(directory+fileName)
        print(sum(ds[('gas', 'cell_mass')]))