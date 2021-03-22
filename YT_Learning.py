# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:22:27 2021

@author: Sherry Wong
"""

import yt
import matplotlib as plt
from yt.units import kpc

import h5py

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# The below code does not work. Is it not an HDF5 file?
# f2 = h5py.File("IsolatedGalaxy/galaxy0030/galaxy0030", 'r')
# print("\nKeys:")
# print(f2.keys())

#create a slice plot
#the vector can be any vector [x, y, z]
#by default the slice is in the middle of the domain
slc = yt.SlicePlot(ds, [0, 0, 1], 'density')
#if you have 2d plots
# slc = yt.plot_2d(ds, 'density')
slc.zoom(100)


#Projection plots of 2 variables?
#https://yt-project.org/doc/visualizing/plots.html#projection-plots
prj = yt.ProjectionPlot(ds, 2, 'temperature', weight_field='density')

# movement around the plotted area
slc.pan((2*kpc, 2*kpc))
# change units and auto-convert
slc.set_axes_unit('Mpc')
# slc.set_axes_unit('Density', 'Msun/pc**3')
# slc.set_center((0.5, 0.53))

###################
#Plot Customization
###################

# # label things
slc.annotate_title("IsolatedGalaxy")
# # you can hide axes/colorbars
# # which you should never do, or Ellen will find you ;)
slc.hide_axes()
# # you can do all sorts of things with your fonts
slc.set_font({'style': 'italic', 'weight': 'bold', 'size': 24})
# # you can change colormaps
slc.set_cmap('density', 'RdBu_r')
# # and use your own custom limits
# slc.set_zlim('density', 10e-30, 10e-25);
# # you can have a log scale (enabled by default)
# use linthresh=num for +/- log scales to adjust near 0.
slc.set_log('density', False)
# # you can set the size of the plot in inches
slc.set_figure_size(10)
# # and resolution
# slc.set_buff_size(1600)
# # you can put stuff on top of your other plots
# # more details at https://yt-project.org/doc/visualizing/callbacks.html#callbacks
slc.annotate_grids()
# # you can also get rid of minor ticks

# # Save plot
# slc.show()
# slc.save("YT_Test_Plots/Customization2")

##################
# 1D Profile Plots
##################

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
my_galaxy = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10*kpc, 3*kpc)

# you can choose any two variables to plot.
# use weight_field to use weighted averages
# syntax ProfilePlot(dataset, x variable, y variable)
plot = yt.ProfilePlot(my_galaxy, 'density', 'temperature')

#you can plot multiple 1D plots 
plot.save("YT_Test_Plots/1DProfile")

################
# 2D Phase Plots
################

# Extra keywords for stats:
# weight_field = none   Total of some field in a bin; default averages values
# fractional = True     Creates probability distribution fcns
# accumulation = True   Cumulative distribution function, i.e. N = sum of 0 to N
phasePlot = yt.PhasePlot(ds, "density", "temperature", ["cell_mass"], weight_field=None)
phasePlot.annotate_title("IsolatedGalaxy Phase Plot weight_field=None")
phasePlot.save("YT_Test_Plots/2DPhasePlot_A")

sphereSelect = ds.sphere('c', 50*kpc)
phasePlot = yt.PhasePlot(sphereSelect, "density", "temperature", ["cell_mass"], fractional=True)
phasePlot.annotate_title("IsolatedGalaxy Phase Plot, 50kpc Sphere Selection")
phasePlot.save("YT_Test_Plots/2DPhasePlot_Sphere")

################
# Particle Plots
################
particlePlot = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y')
particlePlot.annotate_title("IsolatedGalaxy Particle Plot")
particlePlot.save("YT_Test_Plots/2DParticlePlot_A")

# You can add a z axis for a colobar, summed across the box like a projection plot
particlePlot = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y', 'particle_mass')
particlePlot.annotate_title("IsolatedGalaxy Particle Plot + Mass")
particlePlot.save("YT_Test_Plots/2DParticlePlot_B")

################
# Derived Fields
################
def thermal_energy_dens(field, data):
    return (3/2)*data['gas', 'number_density']*data['gas', 'kT']

ds.add_field(('gas', 'thermal_energy_density'),
             units="erg/cm**3", function=thermal_energy_dens,
             sampling_type='cell')

for i in sorted(ds.derived_field_list):
    print(i)
    
ad = ds.all_data
