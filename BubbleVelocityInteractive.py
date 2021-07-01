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
    
    return (posXarray, posYarray, velXarray, velYarray, densityArray, tempArray)
    
    # return {"posXarray" : posXarray, 
    #         "posYarray" : posYarray, 
    #         "velXarray" : velXarray, 
    #         "velYarray" : velYarray, 
    #         "densityArray" : densityArray, 
    #         "tempArray" : tempArray}
    #END OF METHOD
    
#%%
    
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0084"

out = setup(filename)
posXarray = out[0]
posYarray = out[1]
velXarray = out[2]
velYarray = out[3]
densityArray = out[4]
tempArray = out[5]

#%%
ds = yt.load(filename)
ad = ds.all_data()
ymax=float(max(ad['y']).value)
ymin=float(min(ad['y']).value)
x = closestNum(posXarray, 2.32838*pow(10, 22))
ylim = -1.79040182984184e21
slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
#slc.annotate_title("Temperature") #this doesn't work?
#slc.annotate_title("LinePlot")
Xslices = removeDuplicates(posXarray).sort()
# for x in Xslices:
#     yt.linePlot(ds, 'density', [x, 0, 0], [x, ymin, 0], 512)
slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp/0080")

#%%
bounds = {'xmin': 2.1*pow(10, 22), 'xmax': 2.5*pow(10, 22), 'ymin': float(min(ad['y']).value),'ymax': ylim}
x = closestNum(posXarray, 2.32838*pow(10, 22))
#plot many temp plots
filenames = []
for t in range(65, 85):
    filenames.append("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t)) #t is a 2 digit number
    ds = yt.load("/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_00"
                     +str(t))
    #slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    #slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp_linePlot/00"+str(t))
    
    # ad = ds.all_data()
    # dsSelect = ad.include_inside('x', bounds['xmin'], bounds['xmax'])
    # dsSelect = dsSelect.include_inside('y', bounds['ymin'], bounds['ymax'])
    # slc = yt.SlicePlot(ds, 'z', 'temp', data_source=dsSelect, 
    #                center=( np.sum([bounds['xmin'], bounds['xmax']])/2, np.sum([bounds['ymin'], bounds['ymax']])/2, 0))
    # slc.set_width(max([ abs(bounds['xmax']-bounds['xmin']), abs(bounds['ymax']-bounds['ymin']) ]))
    # slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=3,plot_args={'linewidth':0.5,'color':'r'}) 
    # slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp_magStream_slicePlot/00"+str(t))
    
    slc = yt.LinePlot(ds, 'temp', [x, ylim, 0], [x, ymin, 0], 512)
    slc.save("/Users/wongb/Documents/Python_Scripts/YT_Test_Plots/HDF5/temp_linePlots/00"+str(t))

#%%
z=densityArray
lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
# plt.title("Temp (\N{DEGREE SIGN}K)")
plt.clf()
plt.title("Density (g/$cm^3$)")
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values

#%%
#find local extrema
x = closestNum(posXarray, 2.32838*pow(10, 22))
TFtable = (x==posXarray)
tempSliceArray = tempArray[TFtable]
posYSliceArray = posYarray[TFtable]
combo = [(tempSliceArray[0], posYSliceArray[0])] #tie temp points with their Y values
#insert numbers in ascending order
for i in range(1, len(posYSliceArray)-1):
    #print("loop " + str(i))
    #add first term
    # if i == 0: 
    #     combo = np.append(combo, (tempSliceArray[i], posYSliceArray[i]))
    # #add to front
    # elif posYSliceArray[i] < posYSliceArray[0]: 
    #     print("Adding " + (tempSliceArray[i], posYSliceArray[i]) + " at index 0")
    #     combo = np.insert(combo, 0, (tempSliceArray[i], posYSliceArray[i]))
    #add to end
    if combo[len(combo)-1][1] < posYSliceArray[i]: 
        # print("Adding " + str(tempSliceArray[i]) + ", " + str(posYSliceArray[i]) + " at index " + str(i))
        # print("After " + str(combo[len(combo)-1][1]))
        combo.append((tempSliceArray[i], posYSliceArray[i]))
    #add to middle
    else: 
        #print("\t starting inner loop ")
        for x in range(0, len(combo)-1): #sort based on Y value
            # if x == len(combo):
            #     if combo[x][1] > posYSliceArray[i]:
            #         print("problem")
            #     print("Adding (" + str(tempSliceArray[i]) + ", " + str(posYSliceArray[i]) + ") at index " + str(x))
            #     print("After " + str(combo[x][1]))
            #     combo.append((tempSliceArray[i], posYSliceArray[i]))
            #     break
            if combo[x][1] < posYSliceArray[i]:
                # print("Adding (" + str(tempSliceArray[i]) + ", " + str(posYSliceArray[i]) + ") at index " + str(x))
                # print("After " + str(combo[x][1]))
                # print("Before " + str(combo[x+1][1]))
                combo.insert(x+1, (tempSliceArray[i], posYSliceArray[i]))
                break
print(len(combo))
#print(combo)
comboArray = np.asarray(combo)
posYthing = []
numIssues = 0
for i in range(0, len(comboArray)):
    if comboArray[i-1][1] == comboArray[i][1]:
        print(str(i) + ":e " + str(comboArray[i][1]) + ", " + str(comboArray[i][1]) + ", " + str(comboArray[i][1]))
    elif comboArray[i-1][1] > comboArray[i][1]:
        numIssues+=1
        print(str(i) + ":  " + str(comboArray[i][1]) + ", " + str(comboArray[i][1]) + ", " + str(comboArray[i][1]))
print(numIssues)
#%%
extrema = signal.argrelextrema(tempArray, np.greater, order = 2000, axis=0)
print(len(extrema[0]))
print(extrema[0])

extrema = signal.argrelmax(tempArray, np.greater, order = 2000)
print(len(extrema[0]))
print(extrema[0])
