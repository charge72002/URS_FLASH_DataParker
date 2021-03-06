# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:37:53 2021

@author: Sherry Wong
"""

import numpy as np
import yt
from yt.units import kpc
import os

os.chdir("/Users/wongb/Documents/Python Scripts")
filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0020"
#filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"
# filename = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/parkerCRs_hdf5_plt_cnt_0000"

#see https://yt-project.org/doc/analyzing/objects.html
ds = yt.load(filename)

for i in sorted(ds.field_list):
    print(i)
for i in sorted(ds.derived_field_list):
    print(i)
print((('flash', 'cray')) == ('flash', 'cray'));

#####################################
# SlicePlot + Vectors + Streamlines #
#####################################

#Zoom
slc = yt.SlicePlot(ds, 'z', 'density', center=(4.1*kpc, -2.85*kpc, 0), width=(2*kpc, 2*kpc, 0)) #3D!!!
# slc.annotate_velocity(factor=16)
slc.zoom(10);
slc.set_width(2*kpc)
# slc.set_unit('temp', 'K')
slc.annotate_title("parkerCRs_hdf5_plt_cnt_0076 Density + Magnetic Field")
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
# slc.annotate_streamlines('velocity_x','velocity_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
slc.save("YT_Test_Plots/HDF5/StreamMagZoom3");


#No zoom
slc = yt.SlicePlot(ds, 'z', 'density')
# slc.annotate_velocity(factor=16)
slc.annotate_title("parkerCRs_hdf5_plt_cnt_0076 Density + Velocity")
# slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=2,plot_args={'linewidth':0.5,'color':'r'}) 
slc.annotate_streamlines('velocity_x','velocity_y',density=4,plot_args={'linewidth':0.5,'color':'r'}) 
slc.save("YT_Test_Plots/HDF5/StreamVel2")

# slc.hide_axes()
# velMagnitude = np.hypot(ds.r['velx'], ds.r['vely'])
# plot = yt.ProfilePlot(ds, "temperature", ["vely"])
# =============================================================================
# def velMagnitude(field, data):
#     return np.hypot(data['velx'], data['vely'])
# ds.add_field((ds, "velMagnitude"), units="cm/s", function=velMagnitude)
# 

# =============================================================================

##################
# Make phase plots
##################

ad = ds.all_data()
# Use tuples to find individual files; otherwise it will pick a default
plot = yt.PhasePlot(ad, ('gas', "density"), "temperature", ['velocity_magnitude'],
                    fractional = False, accumulation = False)
plot.set_unit('velocity_magnitude', 'cm/s') #vely units weird
plot.annotate_title("Velocity magnitude probability distribution, yes accumulation")
plot.save("YT_Test_Plots/HDF5/TestC")

selection = ds.CutRegion([0 < ('index', 'x') and ('index', 'x') < 2.3*kpc])
plot = yt.SlicePlot(selection, 'z', 'density')
plot.save("YT_Test_Plots/HDF5/VectorSelection");

# Streamlines.streamline(ds, dict['velx'], dict['vely'])
# Streamlines.streamlines.integrate_through_volume()

# print("Saving figures...")
# plt.savefig('Plots/VelocityPlot_m2_c1_0076_stream1.png', dpi=1000);
# #plt.savefig('Plots/VelocityPlot_m2_c1_0076_v3.svg');
# print("Done.")