# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 21:51:17 2021

@author: Sherry Wong; CODE TAKEN FROM VectorPlottingTest.py & ContourPlotting.py
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
import os

os.chdir("/Users/wongb/Documents/Python Scripts")

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
            
    print("Length of posX: " + str(len(posX)))
    posXarray, posYarray = (np.asarray(posX), np.asarray(posY))
    print("Length of posX flattened: " + str(len(posXarray.flatten())))
    posXarray = posXarray.flatten()
    posYarray = posYarray.flatten()
    return posXarray, posYarray
    # END OF METHOD

filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
dict = fileInputTest(filename);
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
velXarray = velXarray.flatten()
velYarray = np.asarray(velY)
velYarray = velYarray.flatten()
print("length of velocity array: " + str(len(velYarray)))
color = np.hypot(velXarray, velYarray)


"""Density Setup""" #from ContourPlotting.py
density = [] #list
#print(dict['density'])
for dens in dict['density']:
    density.append(dens) #shape=(1, 8, 8)
print("length of density list: " + str(len(density)))
densityArray = np.asarray(density)
densityArray = densityArray.flatten()
print("length of new density array: " + str(len(densityArray))) #x64 larger!

#write files for easier access later
with open('m2_c1_0000_data.csv', 'w', newline='') as csvfile:
    fieldnames = ['posx', 'posy', 'velx', 'vely', 'mag', 'density']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    
    for i in range(0, 262143):
        writer.writerow({'posx': posXarray[i], 'posy': posYarray[i], 
                         'velx': velXarray[i], 'vely': velYarray[i], 
                         'mag': color[i], 'density': densityArray[i]})
        #writer.writerow({ posXarray[i], posYarray[i], 
        #                 velXarray[i], velYarray[i], color[i]})

