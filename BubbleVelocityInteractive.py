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

import os
from os import path

import yt
from yt.data_objects.level_sets.api import Clump, find_clumps
yt.toggle_interactivity()

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
    #TODO: if you wanna have fun make this a binary search thing
    for i in smallArray:
        newDev = abs(targetNum-i)
        if(newDev < currentDev):
            currentDev =  newDev
            result = i
    if(result == float("inf")):
        raise Exception("No match found")
    return result
    #END OF METHOD

def setup(fileName):
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
        print("X step: " + str(stepX))
        print("Y step: " + str(stepY))
        
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
    
    dict = fileInputTest(fileName);
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
    
    """"Mag Setup"""
    magX = [] #list
    for thing in dict['magx']:
        magX.append(thing) #shape=(1, 8, 8)
    magXarray = np.asarray(magX)
    magXarray = magXarray.flatten()
    
    magY = [] #list
    for thing in dict['magx']:
        magY.append(thing) #shape=(1, 8, 8)
    magYarray = np.asarray(magY)
    magYarray = magYarray.flatten()

    
    # return (posXarray, posYarray, velXarray, velYarray, densityArray, tempArray, magXarray)
    
    return {"posXarray" : posXarray, 
            "posYarray" : posYarray, 
            "velXarray" : velXarray, 
            "velYarray" : velYarray, 
            "densityArray" : densityArray, 
            "tempArray" : tempArray,
            "magXarray" : magXarray,
            "magYarray" : magYarray}
    #END OF METHOD
    
#%%
### File setup
    
# filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"

out = setup(filename)
posXarray = out["posXarray"]
posYarray = out["posYarray"]
velXarray = out["velXarray"]
velYarray = out["velYarray"]
densityArray = out["densityArray"]
tempArray = out["tempArray"]
magXarray = out["magXarray"]
magYarray = out["magYarray"]

### Plotting setup?
ds = yt.load(filename)
ad = ds.all_data()
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
x = closestNum(posXarray, 2.32838*pow(10, 22))
ylim = -1.79040182984184e21
# slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
# Xslices = removeDuplicates(posXarray).sort()

#%%
### YT Try lineplots a little left or right

dx = posXarray[1]-posXarray[0]
ds = yt.load(filename)
ad = ds.all_data()

x = closestNum(posXarray, 2.32838*pow(10, 22)) - (5*dx) #2.320198554 is good for col 12
for t in range(1, 10):
    line = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    line.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp/0080-"+str(t))
    x = x + dx
    print("("+ str(t) + ", " + str(x) +")")

#%%
### YT annotate lines to visualize where these lineplots are
ylim = -1.79040182984184e21
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])  
slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                   center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
x = closestNum(posXarray, 2.32838*pow(10, 22)) - (5*dx) #2.320198554 is good for col 12
# slc.annotate_grids()
# slc.annotate_cell_edges(line_width = 0.00001, alpha = 0.5)

textx = x + (10**21)
for t in range(1, 10):
    slc.annotate_line((x, -1.8*(10**21), 0), (x, -1.23*(10**22), 0), coord_system="data",  plot_args={"color": "blue"})
    slc.annotate_marker((x, (-2-0.5*t)*(10**21), 0), coord_system="data", plot_args={"color": "red"})
    slc.annotate_text((textx, (-2-0.5*t)*(10**21), 0), str(x), coord_system="data", text_args={"color": "red"})
    x = x + dx
slc.display()

#%%
### YT Plot many temp lineplots over time
ylim = -1.79040182984184e21
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
x = closestNum(posXarray, 2.32838*pow(10, 22)) #A
x = closestNum(posXarray, 2.32019*pow(10, 22)) #B
filenames = []
field = 'temp'
for t in range(65, 85):
    filenames.append("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t)) #t is a 2 digit number
    ds = yt.load("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t))
    
    slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp_linePlotsB/00"+str(t))

#%%
### YT Plot many temp slices over time
ylim = -1.79040182984184e21
saveDirectory = "D:/URS_LargeData/SherryPlots"
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
field = 'temp'

for t in range(65, 85):
    ds = yt.load("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t))
    timeStamp = str(t)
    ad = ds.all_data()
    dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
    dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
        
    slc = yt.SlicePlot(ds, 'z', field,  data_source=dsSelect,
                       center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
    # slc.annotate_velocity(factor=16)
    slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
    slc.annotate_title(timeStamp +" "+ field)
    slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
    slc.set_zlim('temp', 1e2, 1e6)

    if (not path.exists(saveDirectory + "/" + field)):
            os.mkdir(saveDirectory + "/" + field)
            
    slc.save(saveDirectory + "/" + field + "/" + timeStamp)

#%%
### Matplotlib density for interactive window in Spyder
z=densityArray
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.title("Temp (\N{DEGREE SIGN}K)")
plt.clf()
plt.title("Density (g/$cm^3$)")
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

#%%
### Matplotlib temp for interactive window in Spyder 
z=tempArray
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.title("Temp (\N{DEGREE SIGN}K)")
plt.clf()
plt.title("Temp (K)")
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

#%%
### Find local extrema
extrema = signal.argrelextrema(tempArray, np.greater, order = 2000, axis=0)
print(len(extrema[0]))
print(extrema[0])

#%%
### Find local maxima
x = closestNum(posXarray, 2.32019*pow(10, 22)) #B
TFselect = np.logical_and((posXarray == x), (posYarray<0), (posYarray>-1.0*(10**22)))
tempSlice = tempArray[TFselect]
tempSlice = np.array(tempSlice)
# extrema = signal.argrelmax(tempSlice, order = 20)
extrema = signal.find_peaks(tempSlice, height = 2*(10**4))
# extrema = signal.find_peaks(tempSlice, threshold = 2*(10**4))

print("total peaks: " + str(len(extrema[0])))
print("indices:     " + str(extrema[0]))
print("temps:       " + str(tempSlice[extrema[0]]))
print("y positions: " + str(posYarray[extrema[0]]))

### Matplotlib has a mirrored version of the YT plots, but it works
plt.clf()
ySlice = posYarray[TFselect]
# plt.subplot(2, 1, 1)
plt.plot(ySlice, tempSlice)
for y in ySlice[extrema[0]]:
    plt.plot([y, y], [0, 10**6], lw=0.5) #lines [x1, x2] [y1, y2]
plt.yscale('log')
plt.title("local maxima; t=80, x=2.320198554711625e+22")
plt.ylabel("temp(K)")
plt.xlabel("y position(cm)")
# plt.yscale('linear')
plt.show()

#%%
### Maxima: magx inflection points
# ylim = -2.5e21
ylim = -1.79040182984184e21
x = closestNum(posXarray, 2.32019*pow(10, 22)) #B
TFselect = np.logical_and((posXarray == x), (posYarray<ylim), (posYarray>-1.0*(10**22)))
tempSlice = tempArray[TFselect]
ySlice = posYarray[TFselect]
magSlice = magXarray[TFselect]
magSlice = np.array(magSlice)
magDerivative = np.diff(magSlice, axis=0)
magDerivative = np.insert(magDerivative, 0, 0) #to match shape

plt.clf()
# plt.plot(ySlice, magSlice)
# plt.plot(ySlice, magDerivative)
# plt.plot([-(10**22), ylim], [0, 0], lw=0.5) #lines [x1, x2] [y1, y2]
#G = curl([ad['mag_pres'], 0])

# extrema = signal.argrelextrema(magDerivative, comparator=np.greater_equal, order=20)
# extrema = signal.argrelmax(magDerivative, order = 20)
extrema = signal.find_peaks(magDerivative, prominence=(10**-8))
prominences = signal.peak_prominences(magDerivative, extrema[0])
print("total peaks: " + str(len(extrema[0])))
print("indices:     " + str(extrema[0]))
print("magDeriv:    " + str(magDerivative[extrema[0]]))
print("y positions: " + str(posYarray[extrema[0]]))

fig, axs = plt.subplots(2)
fig.suptitle('MagX, MagX Derivative, x='+str(x))
axs[0].plot(ySlice, magSlice)
axs[1].plot(ySlice, magDerivative)
for y in ySlice[extrema[0]]:
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #lines [x1, x2] [y1, y2]
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #prominence lines
    
#%% 
### Experiment with curl?
ds = yt.load(filename)
# ad = ds.all_data()

# slc = yt.SlicePlot(ds, 'z', "magx")
# slc.colorbar()
# slc.display()

def curl(x, y):
    #for our purposes, magX and magY
    xy = np.array([x, y])
    ydx = np.gradient(xy, axis=0)[1]
    xdy = np.gradient(xy, axis=1)[0]
    result = ydx - xdy
    return result

plt.clf
curl = np.asarray(curl(magXarray, magYarray))

ydx = np.gradient([magYarray, posXarray], axis=1)[0]
xdy = np.gradient([magXarray, posYarray], axis=1)[0]
curl = ydx - xdy

# plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
plt.tricontourf(posXarray, posYarray, curl)

#%%
### Experiment with clumps
ds = yt.load(filename)
ylim = -1.79040182984184e21
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
ad = ds.all_data()
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
    
#%%
### Clumps cont'd
c_min = 10 ** np.floor(np.log10(dsSelect['temp']).min())
c_max = 10 ** np.floor(np.log10(dsSelect['temp']).max() + 1)

master_clump = Clump(dsSelect, 'temp')
master_clump.add_validator("min_cells", 20)
## this takes a VERY long time!
find_clumps(master_clump, c_min, c_max, 2.0)
leaf_clumps = master_clump.leaves

prj = yt.ProjectionPlot(ds, "z", 'temp', center="c", width=(20, "kpc"))
prj.annotate_clumps(leaf_clumps)

slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/clump_test_80B")
