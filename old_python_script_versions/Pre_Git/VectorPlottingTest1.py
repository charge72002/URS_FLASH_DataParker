# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:37:53 2021

@author: Sherry Wong
"""

import h5py #installed with "pip install h5py"
import matplotlib.pyplot as plt #Same installation as above
from matplotlib import ticker, cm #for log scale
import matplotlib.colors as colors #for color mapping
# import matplotlib.colors.LogNorm #for logorithmic color mapping
import numpy as np
import sys

def fileInputTest(fileName):
    f2 = h5py.File(fileName, 'r')

    print("\nReading file:")
    print(f2) #general info
    print(f2.attrs)
    #print("\nKeys:")
    #print(f2.keys())

    #print("\nLoading relevant data fields:")
    #print(f2['coordinates'])
    coords = f2['coordinates']

    cray = f2['cray']
    dens = f2['dens']
    magx = f2['magx']
    magy = f2['magy']
    magz = f2['magz']
    pres = f2['pres']
    velx = f2['velx']
    vely = f2['vely']
    velz = f2['velz']

    print() #whitespace
    #print(f2.attrs['Coordinates']).decode("utf-8")

    #for i in f2.attrs["VariableNames"]:
    #    print(i)

    #TODO: Find individual data keys and examine their data types.

    dict = {
        "coordinates" : coords,
        "cray_pressure" : cray,
        "density" : dens,
        "magx" : magx,
        "magy" : magy,
        "magz" : magz,
        "pressure" : pres,
        "velx" : velx,
        "vely" : vely,
        "velz" : velz

    }
    return dict
    #END OF METHOD
    
#filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0020"
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0075"
dict = fileInputTest(filename);


#setup coords
posX = []
posY = []
for coord in dict['coordinates']:
    posX.append(coord[0])
    posY.append(coord[1])
posXarray = np.asarray(posX)
posYarray = np.asarray(posY)
    
#setup velocity
velX = []
velY = []
for x in dict['velx']:
    velX.append(x[0, 0, 0])
    #velX.append(x)
for y in dict['vely']:
    velY.append(y[0, 0, 0])
    #velY.append(y)
velXarray = np.asarray(velX)
velYarray = np.asarray(velY)

plt.quiver(posXarray, posYarray, velXarray, velYarray);
plt.savefig('Plots/VelocityPlot_m2_c1_0075.png', dpi=600); #OK quality
print("Done.")