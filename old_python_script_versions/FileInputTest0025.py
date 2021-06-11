import h5py #installed with "pip install h5py"

filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
f1 = open(filename, 'r')
print(f1) #general info
#print(f1.attrs())

print()

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

filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0025"
dict = fileInputTest(filename)
print(dict)
print(dict['density'])
print(dict['density'][0])

print(". . . . . . . . . . . . . . . . . . . .")
repeat = 'y'
while repeat=='y':
    usrKey = input("Choose a data field. [coordinates/cray_pressure/density/magx/magy/magz/pressure/velx/vely/velz]")
    #print(usrKey)
    usrInt = input("Choose a data point. [0-16383]")
    #print(usrInt)
    print(dict[usrKey])
    #cast as int to avoid "ValueError: Field names only allowed for compound types"
    print(dict[usrKey][int(usrInt)])
    repeat = input("Would you like to see another data field? [y/n]")
print("Goodbye")
#All of the following data sets have 16384 = 2^14 entries b/c 128x128 cells
# print() #<HDF5 dataset "coordinates": shape (16384, 3), type "<f4">
# print(dict['coordinates'])
# print(dict['coordinates'][0])
# for x in [0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16383]:
#     print(dict['coordinates'][x])
# print()
# for x in [100, 1000, 1000, 10000]:
#     print(dict['coordinates'][x])
# print() #<HDF5 dataset "cray": shape (16384, 1, 8, 8), type "<f4">
# print(dict['cray_pressure'])
# print(dict['cray_pressure'][0])
# print() #<HDF5 dataset "dens": shape (16384, 1, 8, 8), type "<f4">
# print(dict['density'])
# print(dict['density'][0])
# print() #<HDF5 dataset "magx": shape (16384, 1, 8, 8), type "<f4">
# print(dict['magx']) #0 for magy and magz
# print(dict['magx'][0])
# print() #<HDF5 dataset "pres": shape (16384, 1, 8, 8), type "<f4">
# print(dict['pressure'])
# print(dict['pressure'][0])
# print() #<HDF5 dataset "vely": shape (16384, 1, 8, 8), type "<f4">
# print(dict['velx']) #0 for velx and velz
# print(dict['velx'][0])
# print(dict['vely']) #0 for velx and velz
# print(dict['vely'][0])
# print(dict['velz']) #0 for velx and velz
# print(dict['velz'][0])

#DO NOT DO THIS, IT WILL PRINT 16384 entries
#for x in dict['density']:
#    print(x)
