import h5py #installed with "pip install h5py"
import matplotlib.pyplot as plt #Same installation as above
import numpy
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
filename = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0025"
# filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0050"
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
    tempXlin = numpy.linspace(coord[0], coord[0] + stepX, 8)
    tempYlin = numpy.linspace(coord[1], coord[1] + stepY, 8)
    tempMeshgrid = numpy.meshgrid(tempXlin, tempYlin)
    posX.append(tempMeshgrid[0])
    posY.append(tempMeshgrid[1])
    #print(coord)
print("Length of posX: " + str(len(posX)))
posXarray, posYarray = (numpy.asarray(posX), numpy.asarray(posY))
print("Length of posX flattened: " + str(len(posXarray.flatten())))

# print(meshgrid)

density = [] #list
print(dict['density'])
for dens in dict['density']:
    density.append(dens)
print("length of density list: " + str(len(density)))
densityArray = numpy.asarray(density)
densityArray = densityArray.flatten()
print("length of new density array: " + str(len(densityArray))) #x64 larger!

fig = plt.figure(figsize=(6, 6))
plt.scatter(posXarray, posYarray,
            linewidths=0, alpha=.7,
            edgecolor='k',
            s = 10,
            c=densityArray)
#plt.show()
#^ from https://stackoverflow.com/questions/59232073/scatter-plot-with-3-variables-in-matplotlib
print("\nPlotting...")
plt.savefig('Plots/CoordDensPlot_m1.5_c1_0025.png'); #bad quality
plt.savefig('Plots/CoordDensPlot_m1.5_c1_0025.svg'); #much better quality
print("\nPlotting finished.")
