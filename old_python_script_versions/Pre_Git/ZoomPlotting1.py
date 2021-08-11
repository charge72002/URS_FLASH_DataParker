# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 23:16:05 2021

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
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0076"
dict = fileInputTest(filename);
        
"""Sets up coordinates for 16384 * 64 data points. Returns posXarray, posYarray"""
def bigCoordinateSetup():
    #Make lists of individual data columns
    posX = [] #list
    posY = [] #list
    
    #xOffset = int()
    stepX = dict['coordinates'][1][0] - dict['coordinates'][0][0]
    stepY = dict['coordinates'][2][1] - dict['coordinates'][1][1] #bug was here
    
    for coord in dict['coordinates']:
        tempXlin = np.linspace(coord[0], coord[0] + stepX, 8)
        tempYlin = np.linspace(coord[1], coord[1] + stepY, 8)
        tempMeshgrid = np.meshgrid(tempXlin, tempYlin)
        posX.append(tempMeshgrid[0])
        posY.append(tempMeshgrid[1]) #bug is not here
        #print(coord)
            
    # print("Length of posX: " + str(len(posX)))
    posXarray, posYarray = (np.asarray(posX), np.asarray(posY))
    # print("Length of posX flattened: " + str(len(posXarray.flatten())))
    posXarray = posXarray.flatten()
    posYarray = posYarray.flatten()
    return posXarray, posYarray
    # END OF METHOD
posXarray, posYarray = bigCoordinateSetup()

"""Velocity Setup"""
velX = []
velY = []
for x in dict['velx']:
    #velX.append(x[0, 0, 0])
    velX.append(x)
for y in dict['vely']:
    #velY.append(y[0, 0, 0])
    velY.append(y)
velXarray = np.asarray(velX)
velYarray = np.asarray(velY)

velXarray = velXarray.flatten()
velYarray = velYarray.flatten()

"""Density Setup"""
density = [] #list
#print(dict['density'])
for dens in dict['density']:
    density.append(dens) #shape=(1, 8, 8)
densityArray = np.asarray(density)
densityArray = densityArray.flatten()

# """Adjust boundaries"""
TFtable = np.logical_and((pow(10, 22) < posXarray), (posYarray < 0))

posXarray = posXarray[TFtable]
posYarray = posYarray[TFtable]
densityArray = densityArray[TFtable]
velXarray = velXarray[TFtable]
velYarray = velYarray[TFtable]


print("Plotting...")


"""Plot Density"""
lev = np.logspace(np.log10(densityArray.min()), np.log10(densityArray.max()), num=100)
plt.tricontourf(posXarray, posYarray, densityArray, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

upperTick = round( np.log10(densityArray.max()) )
bottomTick = round( np.log10(densityArray.min()) )
t=np.logspace(upperTick, bottomTick, round( int(upperTick-bottomTick)/2) )
plt.colorbar(extend='both', orientation='vertical', ticks=t)

"""Plot Velocity"""
# color = np.hypot(velXarray, velYarray)
# plt.quiver(posXarray, posYarray, velXarray, velYarray, velYarray, norm=matplotlib.colors.LogNorm(),
#            )
#            #width=0.001, minlength = 0.001, headwidth=0.25, headlength=0.5)
# plt.colorbar(extend='both', orientation='vertical')

#plt.scatter(posXarray, posYarray, color='0.5', s=1)

print("Saving figures...")
# plt.savefig('Plots/DensityPlot_m2_c1_0076_res5.png', dpi=1000);
plt.savefig('Plots/DensityPlot_m2_c1_0076_zoom4.png', dpi=1000);
#plt.savefig('Plots/VelocityPlot_m2_c1_0076_v3.svg');
print("Done.")