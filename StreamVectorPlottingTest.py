# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:37:53 2021

@author: Sherry Wong
"""

import h5py #installed with "pip install h5py"
import matplotlib
import matplotlib.pyplot as plt #Same installation as above
from matplotlib import ticker, cm #for log scale
import matplotlib.colors as colors #for color mapping
# import matplotlib.colors.LogNorm #for logorithmic color mapping
import numpy as np
import sys
import csv #https://docs.python.org/3/library/csv.html#csv.Dialect
from yt.visualization.api import Streamlines

import yt
from yt.units import kpc

#filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0020"
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"

#see https://yt-project.org/doc/analyzing/objects.html
ds = yt.load(filename)

#Zoom
slc = yt.SlicePlot(ds, 'z', 'density', center=(4*kpc, -2.5*kpc, 0), width=(2*kpc, 2*kpc, 0)) #3D!!!
slc.annotate_velocity(factor=16)
slc.zoom(10);
slc.set_width(2*kpc)
# slc.set_unit('temp', 'K')
slc.annotate_title("parkerCRs_hdf5_plt_cnt_0076 Density + Velocity")
slc.save("YT_Test_Plots/HDF5/VectorZoom2");

#No zoom
slc = yt.SlicePlot(ds, 'z', 'density')
slc.annotate_velocity(factor=16)
slc.annotate_title("parkerCRs_hdf5_plt_cnt_0076 Density + Velocity")
slc.save("YT_Test_Plots/HDF5/Vector2");

# slc.hide_axes()
# velMagnitude = np.hypot(ds.r['velx'], ds.r['vely'])
# plot = yt.ProfilePlot(ds, "temperature", ["vely"])
# =============================================================================
# def velMagnitude(field, data):
#     return np.hypot(data['velx'], data['vely'])
# ds.add_field((ds, "velMagnitude"), units="cm/s", function=velMagnitude)
# 
# for i in sorted(ds.field_list):
#     print(i)
# for i in sorted(ds.derived_field_list):
#     print(i)
# =============================================================================

##################
# Make phase plots
##################

ad = ds.all_data()
plot = yt.PhasePlot(ad, "density", "temperature", ['velocity_magnitude'],
                    fractional = True, accumulation = True)
plot.set_unit('velocity_magnitude', 'cm/s') #vely units weird
plot.annotate_title("Velocity magnitude probability distribution, yes accumulation")
plot.save("YT_Test_Plots/HDF5/magB")

# Streamlines.streamline(ds, dict['velx'], dict['vely'])
# Streamlines.streamlines.integrate_through_volume()

# print("Saving figures...")
# plt.savefig('Plots/VelocityPlot_m2_c1_0076_stream1.png', dpi=1000);
# #plt.savefig('Plots/VelocityPlot_m2_c1_0076_v3.svg');
# print("Done.")