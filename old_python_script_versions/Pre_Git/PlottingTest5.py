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

#filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
# filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0110"
# filename = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/parkerCRs_hdf5_plt_cnt_0050"
dict = fileInputTest(filename)

print(". . . . . . . . . . . . . . . . . . . .")
#All of the following data sets have 16384 = 2^14 entries b/c 128x128 cells

#plt.plot([1, 2, 3, 4], [1, 2, 3, 4])
#print("\nPlotting finished.")

#Make lists of individual data columns
posX = [] #list
posY = [] #list

#xOffset = int()
stepX = dict['coordinates'][1][0] - dict['coordinates'][0][0]
stepY = dict['coordinates'][1][1] - dict['coordinates'][0][1]

for coord in dict['coordinates']:
    tempXlin = np.linspace(coord[0], coord[0] + stepX, 8)
    tempYlin = np.linspace(coord[1], coord[1] + stepY, 8)
    tempMeshgrid = np.meshgrid(tempXlin, tempYlin)
    posX.append(tempMeshgrid[0])
    posY.append(tempMeshgrid[1])
    #print(coord)
print("Length of posX: " + str(len(posX)))
posXarray, posYarray = (np.asarray(posX), np.asarray(posY))
print("Length of posX flattened: " + str(len(posXarray.flatten())))
posXarray, posYarray = (posXarray.flatten(), posYarray.flatten())
# print(meshgrid)

density = [] #list
#print(dict['density'])
for dens in dict['density']:
    density.append(dens)
print("length of density list: " + str(len(density)))
densityArray = np.asarray(density)
densityArray = densityArray.flatten()
print("length of new density array: " + str(len(densityArray))) #x64 larger!

# fig = plt.figure(figsize=(6, 6))
# pcm = plt.pcolor(posXarray, posYarray, densityArray,
#            norm=colors.LogNorm(vmin=densityArray.min(), vmax=densityArray.max()),
#            cmap='PuBu_r')
# fig.colorbar(pcm, extend='max')
# plt.title("Density")
# plt.scatter(posXarray, posYarray,
#              linewidths=0, alpha=.7,
#              edgecolor='k',
#              s = 2,
#              cmap=cm.PuBu_r,
#              c=densityArray)

# print("reshaping...")
# rowcol = int(np.sqrt(len(densityArray))) #1D to 2D helper
# X, Y = np.reshape(posXarray, (rowcol, rowcol)), np.reshape(posYarray, (rowcol, rowcol))
# Z = np.reshape(densityArray, (rowcol, rowcol))
# print(X.shape)
# print(Y.shape)
# print(Z.shape)
# print(Z.min())
# print(Z.max())
# print("reshaping done!")

# print("\nPlotting...")

# pcm = plt.pcolormesh(X, Y, Z, #min has be increased to 1E-41
#                    norm=colors.LogNorm(vmin=Z.min()+0.000000000000000000000000000000000000001, vmax=Z.max()),)
# plt.colorbar(pcm, extend='max')
# plt.title('Density')

# fig, ax = plt.subplots(2, 1)

# pcm = ax[0].pcolormesh(X, Y, Z,
#                    norm=colors.LogNorm(vmin=Z.min()+0.0000001, vmax=Z.max()),
#                    cmap='PuBu_r')
# fig.colorbar(pcm, ax=ax[0], extend='max')

# pcm = ax[1].pcolormesh(X, Y, Z, cmap='PuBu_r')
# fig.colorbar(pcm, ax=ax[1], extend='max')

# plt.show()


#^ from https://matplotlib.org/3.2.1/tutorials/colors/colormapnorms.html

#trying different plotting methods:
#XY = np.asarray([posXarray, posYarray])
plt.tricontourf(posXarray, posYarray, densityArray, levels = 50,
                norm=colors.LogNorm(vmin=densityArray.min(), vmax=densityArray.max()))
                #cmap=colors.Colormap('PuBu_r', N=50))
plt.colorbar()
#plt.pcolormesh(XY, densityArray)
#plt.show()
# ^ from https://stackoverflow.com/questions/59232073/scatter-plot-with-3-variables-in-matplotlib

#plt.savefig('Plots/CoordDensPlot_m1.5_c1_0025.png'); #bad quality
plt.savefig('Plots/CoordDensPlot_m1.5_c1_0000.png', dpi=600); #OK quality
#plt.savefig('Plots/CoordDensPlot_m2_c1_0111.svg'); #much better quality, but uses a LOT of space
print("\nPlotting finished.")
