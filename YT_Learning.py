# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:22:27 2021

@author: Sherry Wong
"""

import yt
import matplotlib as plt
from yt.units import kpc

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

#create a slice plot
#the vector can be any vector [x, y, z]
slc = yt.SlicePlot(ds, [1, 1, 1], 'density')
#if you have 2d plots
# slc = yt.plot_2d(ds, 'density')
slc.zoom(100)


#Projection plots of 2 variables?
#https://yt-project.org/doc/visualizing/plots.html#projection-plots
prj = yt.ProjectionPlot(ds, 2, 'temperature', weight_field='density')



#Save plot
slc.show()
slc.save("YT_Test_Plots/Slice4")
