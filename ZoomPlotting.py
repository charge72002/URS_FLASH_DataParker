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

""""Temp Setup"""
temp = [] #list
for thing in dict['temp']:
    temp.append(thing) #shape=(1, 8, 8)
tempArray = np.asarray(temp)
tempArray = tempArray.flatten()

"""Adjust boundaries"""
TFtable = np.logical_and((pow(10, 22) < posXarray), (posYarray < 0))

# mag = np.hypot(velXarray, velYarray)
# TFtable = np.logical_and(5*pow(10, 6) < mag))

conversion = 3.086e21 #kpc to cm
#TFtable = np.logical_and((-1.5*conversion<posXarray),(posXarray<0.5*conversion), (posYarray < 0))

posXarray = posXarray[TFtable]
posYarray = posYarray[TFtable]
densityArray = densityArray[TFtable]
tempArray = tempArray[TFtable]
velXarray = velXarray[TFtable]
velYarray = velYarray[TFtable]

#write files for easier access later
# with open('m2_c1_0076_data3.csv', 'w', newline='') as csvfile:
#     fieldnames = ['posx', 'posy', 'velx', 'vely', 'mag']
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#     writer.writeheader()
    
#     color = np.hypot(velXarray, velYarray)

#     for i in range(0, 262143):
#         writer.writerow({'posx': posXarray[i], 'posy': posYarray[i], 
#                           'velx': velXarray[i], 'vely': velYarray[i], 'mag': color[i]})
#         #writer.writerow({ posXarray[i], posYarray[i], 
#         #                 velXarray[i], velYarray[i], color[i]})

print("Plotting...")


"""Plot Density"""
z=densityArray
# plt.title("Temp (\N{DEGREE SIGN}K)")
plt.title("Density (g/$cm^3$)")

lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
#plt.tricontourf(posXarray, posYarray, z, locator=ticker.LogLocator())
                
upperTick = round( np.log10(z.max()) )
bottomTick = round( np.log10(z.min()) )

numTicks = upperTick-bottomTick
if numTicks > 5:
    numTicks= 6+(upperTick-bottomTick)%5
else:
    numTicks+=1


# t = np.logspace( upperTick, bottomTick, round( int( (upperTick-bottomTick)) )+1 )
t = np.logspace( upperTick, bottomTick, int(numTicks) )
t = np.insert(t, 0, z.min())
t = np.insert(t, len(t), z.max())
print("min, max = " + str(z.min()) + ", " + str(z.max()))
# t=np.asarray(t)
# t=t.flatten()
plt.colorbar(extend='both', orientation='vertical', ticks=t, format="%1.2g")
# see https://docs.python.org/2/library/stdtypes.html#string-formatting ^

"""Plot Velocity"""
# color = np.hypot(velXarray, velYarray)
# plt.quiver(posXarray, posYarray, velXarray, velYarray, color, norm=matplotlib.colors.LogNorm(),
#             cmap = plt.get_cmap('copper'))
#             #width=0.001, minlength = 0.001, headwidth=0.25, headlength=0.5)
# plt.colorbar(extend='both', orientation='vertical')

#plt.scatter(posXarray, posYarray, color='0.5', s=1)

plt.xlabel("x (kpc)")
plt.ylabel("y (kpc)")

print("Saving figures...")
# plt.savefig('Plots/DensityPlot_m2_c1_0076_zoom5.png', dpi=1000);
# plt.savefig('Plots/TempPlot_m2_c1_0076_res5.png', dpi=1000);
#plt.savefig('Plots/VelocityPlot_m2_c1_0076_v3.svg');
print("Done.")