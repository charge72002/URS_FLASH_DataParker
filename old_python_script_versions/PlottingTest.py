import h5py #installed with "pip install h5py"
import matplotlib.pyplot as plt #Same installation as above
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

    #print(f2['cray'])
    cray = f2['cray']

    #print(f2['dens'])
    dens = f2['dens']

    #print(f2['magx'])
    magx = f2['magx']
    #print(f2['magy'])
    magy = f2['magy']
    #print(f2['magz'])
    magz = f2['magz']

    #print(f2['pres'])
    pres = f2['pres']

    #print(f2['velx'])
    velx = f2['velx']
    #print(f2['vely'])
    vely = f2['vely']
    #print(f2['velz'])
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

print("\nPrinting dictionaries . . . . . . . . . . . . . . . . . . . .")

#Working Directory
# pwd = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
# fileName = input("Enter a file name.")
# filePath = pwd + fileName
# print(filePath)

#filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
dict = fileInputTest(filename)
print(dict)
print(dict['density'])
print(dict['density'][0])

print(". . . . . . . . . . . . . . . . . . . .")
#All of the following data sets have 16384 = 2^14 entries b/c 128x128 cells
# =============================================================================
# print() #<HDF5 dataset "coordinates": shape (16384, 3), type "<f4">
# print(dict['coordinates'])
# print(dict['coordinates'][0])
# print() #<HDF5 dataset "dens": shape (16384, 1, 8, 8), type "<f4">
# print(dict['density'])
# print(dict['density'][0])
# =============================================================================

#plt.plot([1, 2, 3, 4], [1, 2, 3, 4])
#print("\nPlotting finished.")

#Make lists of individual data columns
posX = [] #list
posY = [] #list
for coord in dict['coordinates']:
    posX.append(coord[0])
    #xOffset = int()
    #posX.append([[0, 1, 2, 3, 4, 5, 6, 7], [0, 1, 2, 3, 4, 5, 6, 7]])
    posY.append(coord[1])
    print(coord)
density = [] #list
#print(dict['density'])
for dens in dict['density']:
    density.append(dens[0][4][4])
#print(len(density))

#recreate the z-order curve:
# =============================================================================
# plt.plot(posX, posY, linewidth = 0.1)
# plt.savefig('CoordPlot1.png'); #bad quality
# plt.savefig('CoordPlot1.svg'); #much better quality
# print("\nPlotting finished.")
# =============================================================================

#plot coordinates v density
#plt.plot(posX, posY, density, s=500, cmap='gray')
#^ from https://stackoverflow.com/questions/8202605/matplotlib-scatterplot-colour-as-a-function-of-a-third-variable
fig = plt.figure(figsize=(6, 6))
plt.scatter(posX, posY,
           linewidths=0, alpha=.7,
           edgecolor='k',
           s = 5,
           c=density)
#plt.show()
#^ from https://stackoverflow.com/questions/59232073/scatter-plot-with-3-variables-in-matplotlib
plt.savefig('CoordDensPlot2.png'); #bad quality
plt.savefig('CoordDensPlot2.svg'); #much better quality
print("\nPlotting finished.")
