#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 10:55:40 2021

@author: Sherry Wong

Clumps:
https://yt-project.org/doc/analyzing/domain_analysis/clump_finding.html
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
#%%
import numpy as np

import yt
from yt.data_objects.level_sets.api import *
from yt.data_objects.level_sets.api import Clump
yt.toggle_interactivity()

import matplotlib.pyplot as plt
from matplotlib import ticker #for log scale

import sys
sys.path.insert(0, pwd)
# Sherry packages
# import estimate_bubbleVelocity as est
import hdf5_parser
#%%
ds = yt.load(filename)
ad = ds.all_data()
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
ylim = -1.79040182984184e21

bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])

#%%
c_min = 10 ** np.floor(np.log10(dsSelect['density']).min())
c_max = 10 ** np.floor(np.log10(dsSelect['density']).max() + 1)

master_clump = Clump(dsSelect, 'density')

#satisfy this condition to be your own clump
#pre-made validators are minimum cells and gravitational boundedness
master_clump.add_validator("min_cells", 20)

#add a custom validator
def _minimum_temp(clump, min_temp):
    return False not in (clump['temp'] >= min_temp)
add_validator("min_temp", _minimum_temp)
master_clump.add_validator("min_temp", (10**4))

#if you get keyerror 1350, rerun/setup the code from the top.
find_clumps(master_clump, c_min, c_max, 20.0)
leaf_clumps = master_clump.leaves

#%%
# prj = yt.ProjectionPlot(dsSelect, "z", 'temp')
# prj.annotate_clumps(leaf_clumps)
slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                       center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
# slc.annotate_velocity(factor=16)
# slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
slc.annotate_title("Temp clumps, temp>10^4K, min_cells = 20")
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
slc.annotate_clumps(leaf_clumps)
slc.save(pwd + "/Plots/clumps")

#%%
out = hdf5_parser.setup(filename)
z=out['densityArray']
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.title("Temp (\N{DEGREE SIGN}K)")
# plt.close()
plt.clf()
fig, ax = plt.subplots(1)
ax[0].tricontourf(out['posXarray'], out['posYarray'], z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
ax[0].title("Density (g/$cm^3$)")
cbar = plt.colorbar(extend='both', orientation='vertical')
cbar.ax.set_yticklabels(lev)