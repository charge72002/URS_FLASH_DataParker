# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 19:21:31 2021

@author: Sherry Wong
"""
import h5py
import numpy as np
import scipy
from scipy import signal

import matplotlib
import matplotlib.pyplot as plt #Same installation as above
from matplotlib import ticker, cm #for log scale
import matplotlib.colors as colors #for color mapping

import yt

filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0075"

def setup:
    def fileInputTest(fileName):
        f2 = h5py.File(fileName, 'r')
    
        print("\nReading file:")
        print(f2) #general info
        print(f2.attrs)
        print("\nKeys:")
        print(f2.keys())
    
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
        temp = f2['temp']
    
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
            "velz" : velz,
            "temp" : temp
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
                
        # print("Length of posX: " + str(len(posX)))
        posXarray, posYarray = (np.asarray(posX), np.asarray(posY))
        # print("Length of posX flattened: " + str(len(posXarray.flatten())))
        posXarray = posXarray.flatten()
        posYarray = posYarray.flatten()
        return posXarray, posYarray
        # END OF METHOD
    
    def removeDuplicates(arrayIN):
        listOUT = []
        for i in arrayIN:
            if(listOUT.count(i)==0):
                listOUT.append(i)        
        arrayOUT = np.asarray(listOUT)
        return arrayOUT
        #END OF METHOD
    
    #Finds the number in the array closest to the target number.
    def closestNum(arrayIN, targetNum):
        smallArray = removeDuplicates(arrayIN)
        currentDev = float("inf")
        result = float("inf")
        for i in smallArray:
            newDev = abs(targetNum-i)
            if(newDev < currentDev):
                currentDev =  newDev
                result = i
        if(result == float("inf")):
            raise Exception("No match found")
        return result
        #END OF METHOD
    
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
    
    """"Temp Setup"""
    temp = [] #list
    for thing in dict['temp']:
        temp.append(thing) #shape=(1, 8, 8)
    tempArray = np.asarray(temp)
    tempArray = tempArray.flatten()
    #END OF METHOD

ds = yt.load(filename)
ad = ds.all_data()
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
x = 2.32838*pow(10, 22)
ylim = -1.79040182984184e21
slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
slc.annotate_title("Temperature")
#slc.annotate_title("LinePlot")
Xslices = removeDuplicates(posXarray).sort()
# for x in Xslices:
#     yt.linePlot(ds, 'density', [x, 0, 0], [x, ymin, 0], 512)
slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp/0075")

#plot many temp plots
filenames = []
for t in range(65, 85):
    filenames.append("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t)) #t is a 2 digit number
    ds = yt.load("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t))
    slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp/00"+str(t))

    
z=densityArray
# plt.title("Temp (\N{DEGREE SIGN}K)")
plt.title("Density (g/$cm^3$)")

lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

#find local extrema
TFtable = (x==posXarray)
tempSliceArray = tempArray[TFtable]
tempX = np.array([tempArray, posYarray])
# extrema = signal.argrelextrema(tempArray, np.greater, order = 2000, axis=0)
# print(len(extrema[0]))
# print(extrema[0])
