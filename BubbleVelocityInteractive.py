# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 19:21:31 2021

@author: Sherry Wong
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
# import h5py
import numpy as np
# import scipy
from scipy import signal

# import matplotlib
import matplotlib.pyplot as plt #Same installation as above
from matplotlib import ticker, cm #for log scale
import matplotlib.colors as colors #for color mapping

import os
from os import path
import csv

import yt
from yt.data_objects.level_sets.api import Clump, find_clumps
yt.toggle_interactivity()

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import hdf5_parser

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

## returns the indices right before a sign flip
## i.e. findSignFlips([-3, -1, 1, 3, -2]) returns [1, 3]
def findSignFlips(arrayIN):
    result = []
    prevNum = arrayIN[0]
    for index in range(0, len(arrayIN)):
       if prevNum < 0 and arrayIN[index] > 0: result.append(index-1)
       elif prevNum > 0 and arrayIN[index] < 0: result.append(index-1)
       prevNum = arrayIN[index]
    return np.array(result)
#%%
### File setup
out = hdf5_parser.setup(filename)
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
dx = posXarray[1]-posXarray[0]
# slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
# Xslices = removeDuplicates(posXarray).sort()

#%%
### YT Try lineplots a little left or right
ds = yt.load(filename)
ad = ds.all_data()

x = closestNum(posXarray, 2.32838*pow(10, 22)) - (5*dx) #2.320198554 is good for col 12
for t in range(1, 10):
    line = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    line.save(pwd + "/YT_Test_Plots/HDF5/temp/0080-"+str(t))
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
slc.save(pwd + "/bubble_velocity/line_annotations")
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
    filenames.append(filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)) #t is a 2 digit number
    ds = yt.load(filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t))
    
    slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp_linePlotsB/00"+str(t))

#%%
### YT Plot many temp slices over time
ylim = -1.79040182984184e21
saveDirectory = "D:/URS_LargeData/SherryPlots"
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
field = 'temp'

for t in range(65, 85):
    ds = yt.load(filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t))
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
# plt.close()
plt.clf()
plt.title("Density (g/$cm^3$)")
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

#%%
### Matplotlib temp for interactive window in Spyder 
z=tempArray
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.title("Temp (\N{DEGREE SIGN}K)")
# plt.close()
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
ySlice = out['posYarray'][TFselect]
tempSlice = tempArray[TFselect]
tempSlice = np.array(tempSlice)
# extrema = signal.argrelmax(tempSlice, order = 20)
extrema = signal.find_peaks(tempSlice, height = 2*(10**4))
# extrema = signal.find_peaks(tempSlice, threshold = 2*(10**4))

print("total peaks: " + str(len(extrema[0])))
print("indices:     " + str(extrema[0]))
print("temps:       " + str(tempSlice[extrema[0]]))
print("y positions: " + str(ySlice[extrema[0]]))

### Matplotlib has a mirrored version of the YT plots, but it works
# plt.close()
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

t=80
filenameB = filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)
out = hdf5_parser.setup(filenameB)
    
# ylim = -2.5e21
ylim = -1.79040182984184e21
x = closestNum(out['posXarray'], 2.32019*pow(10, 22)) #B
TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>-1.0*(10**22)))

tempSlice = out['tempArray'][TFselect]
ySlice = out['posYarray'][TFselect]
magSlice = np.array( out['magXarray'][TFselect] )
magDerivative = np.diff(magSlice, axis=0)
magDerivative = np.insert(magDerivative, 0, 0) #to match shape

# plt.close()
plt.clf()
# plt.plot(ySlice, magSlice)
# plt.plot(ySlice, magDerivative)
# plt.plot([-(10**22), ylim], [0, 0], lw=0.5) #lines [x1, x2] [y1, y2]
#G = curl([ad['mag_pres'], 0])

# extrema = signal.argrelextrema(magDerivative, comparator=np.greater_equal, order=20)
# extrema = signal.argrelmax(magDerivative, order = 20)
extrema = signal.find_peaks(magDerivative, prominence=(5*(10**-9))) #height=6*(10**-9)
prominences = signal.peak_prominences(magDerivative, extrema[0])
print("not adjusted")
print("total peaks: " + str(len(extrema[0])))
print("indices:     " + str(extrema[0]))
print("magDeriv:    " + str(magDerivative[extrema[0]]))
print("y positions: " + str(ySlice[extrema[0]]))
print("prominences: " + str(prominences))

## subtract the initial condition
zeroFile = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
zeroOut = hdf5_parser.setup(zeroFile)

magXarrayAdj = out['magXarray'] - zeroOut['magXarray']
magSliceAdj = np.array( magXarrayAdj[TFselect] )
magDerivativeAdj = np.diff(magSlice, axis=0)
magDerivativeAdj = np.insert(magDerivativeAdj, 0, 0) #to match shape
#%%
extremaAdj = signal.find_peaks(magDerivativeAdj, prominence=(5*(10**-9)), width=20) #height=6*(10**-9)
prominencesAdj = signal.peak_prominences(magDerivativeAdj, extrema[0])
print("adjusted")
print("total peaks: " + str(len(extremaAdj[0])))
print("indices:     " + str(extremaAdj[0]))
print("magDeriv:    " + str(magDerivativeAdj[extremaAdj[0]]))
print("y positions: " + str(ySlice[extremaAdj[0]]))
print("prominences: " + str(prominencesAdj))

zerosTF = np.logical_and(magSliceAdj<10**-10, -10**-10<magSliceAdj)
magXarrayZeros = magSliceAdj[zerosTF]
ySliceZeros = ySlice[zerosTF]
print(ySliceZeros)
print(magXarrayZeros.shape)

fig, axs = plt.subplots(4)

axs[0].plot(ySlice, magSlice)
axs[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
for y in ySliceZeros:
    axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #lines [x1, x2] [y1, y2]
    axs[0].plot([y, y], [np.min(magSlice), np.max(magSlice)], lw=1) #prominence lines
axs[1].plot(ySlice, magDerivative)
for y in ySlice[extrema[0]]:
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #lines [x1, x2] [y1, y2]
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #prominence lines
    
axs[2].plot(ySlice, magSliceAdj)
axs[2].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
for y in ySliceZeros:
    axs[2].plot([y, y], [np.min(magSliceAdj), np.max(magSliceAdj)], lw=1) #lines [x1, x2] [y1, y2]
    axs[2].plot([y, y], [np.min(magSliceAdj), np.max(magSliceAdj)], lw=1) #prominence lines
axs[3].plot(ySlice, magDerivativeAdj)
for y in ySlice[extremaAdj[0]]:
    axs[3].plot([y, y], [0, np.max(magDerivativeAdj)], lw=1) #lines [x1, x2] [y1, y2]
    axs[3].plot([y, y], [0, np.max(magDerivativeAdj)], lw=1) #prominence lines
    
fig.savefig(pwd + "/bubble_velocity/AdjustedB_t" + str(t))

#%%
fig, axs = plt.subplots(2)
fig.suptitle('MagX, MagX Derivative, x='+str(x))
axs[0].plot(ySlice, magSlice)
axs[0].plot([max(ySlice), min(ySlice)], [0, 0], lw=1)
axs[1].plot(ySlice, magDerivative)
for y in ySlice[extrema[0]]:
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #lines [x1, x2] [y1, y2]
    axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #prominence lines
fig.savefig("/Users/bwong/URS_FLASH_DataParker/bubble_velocity/col12_t" + str(t))

#%%
# matplotlib plot position onto full, large slice
fig, axs = plt.subplots(1)
axs.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
for y in ySlice[extrema[0]]:
    axs.plot([x, x],[bounds['ymin'], bounds['ymax']], lw=1, color='red')
    axs.plot([bounds['xmin'], bounds['xmax']], [y, y], lw=1, color='red')

#%%
# YT plot position onto slice
# select col 12
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
ad = ds.all_data()
dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
z=tempArray
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
slc = yt.SlicePlot(ds, 'z', 'temp',  data_source=dsSelect,
                       center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
# slc.annotate_title(timeStamp +" "+ field)
slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
slc.set_zlim('temp', 1e2, 1e6)
# ann_slc = plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev)
# plt.title("Temp (\N{DEGREE SIGN}K)")
for y in ySlice[extrema[0]]:
    slc.annotate_line((x, bounds['ymin'], 0), (x, bounds['ymax'], 0), coord_system="data",  plot_args={"color": "red", "linewidth": 1})
    slc.annotate_line((bounds['xmin'], y, 0), (bounds['xmax'], y, 0), coord_system="data",  plot_args={"color": "red",  "linewidth": 1})
slc.save(pwd + "/bubble_velocity/col12_t80")
slc.display()
#%%
### Iterate peak detection through time
total_peaks = []
Ypos = []
prominences = []
zeroFile = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
zeroOut = hdf5_parser.setup(zeroFile)
x = closestNum(out['posXarray'], 2.32019*pow(10, 22))
for t in range(65, 85):
    filename = filedirectory + "/parkerCRs_hdf5_plt_cnt_00" + str(t)
    out = hdf5_parser.setup(filename)
    
    #set up fields
    TFselect = np.logical_and((out['posXarray'] == x), (out['posYarray']<ylim), (out['posYarray']>-1.0*(10**22)))
    magX = out['magXarray'] - zeroOut['magXarray']
    ySlice = out['posYarray'][TFselect]
    magSlice = magX[TFselect]
    magSlice = np.array(magSlice)
    magDerivative = np.diff(magSlice, axis=0)
    magDerivative = np.insert(magDerivative, 0, 0) #to match shape
    
    extrema = signal.find_peaks(magDerivative, prominence=(5*(10**-9)))
    total_peaks.append(len(extrema[0]))
    Ypos.append(ySlice[extrema[0]])
    prominences.append(signal.peak_prominences(magDerivative, extrema[0])[0])
    
    fig, axs = plt.subplots(2)
    fig.suptitle('MagX, MagX Derivative, t='+ str(t) +', x='+str(x))
    axs[0].plot(ySlice, magSlice)
    axs[0].plot([np.max(ySlice), np.min(ySlice)], [0, 0], lw=1) #0 line
    axs[1].plot(ySlice, magDerivative)
    for y in ySlice[extrema[0]]:    
        axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #lines [x1, x2] [y1, y2]
        axs[1].plot([y, y], [0, np.max(magDerivative)], lw=1) #prominence lines
    fig.savefig(pwd + "/bubble_velocity/col12_t" + str(t))
    plt.close()

#%%
### save important data as .csv
with open(pwd + '/bubble_velocity/col12.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', '#peaks', 'Ypos', 'prominences']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(0, len(total_peaks)):
        writer.writerow({'t': i+65, 
                         '#peaks': total_peaks[i],
                         'Ypos': Ypos[i], 
                         'prominences': prominences[i]})
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

# plt.close()
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
