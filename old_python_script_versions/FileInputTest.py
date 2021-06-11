import h5py #installed with "pip install h5py"

filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
f1 = open(filename, 'r')
print(f1) #general info
#print(f1.attrs())

print()

def fileInputTest(fileName):
    f2 = h5py.File(fileName, 'r')
    print(f2) #general info
    print(f2.attrs)
    print("\nKeys:")
    print(f2.keys())

    print("\nLoading relevant data fields:")
    print(f2['coordinates'])
    coords = f2['coordinates']

    print(f2['dens'])
    dens = f2['dens']

    print(f2['magx'])
    magx = f2['magx']
    print(f2['magy'])
    magy = f2['magy']
    print(f2['magz'])
    magz = f2['magz']

    print(f2['pres'])
    pres = f2['pres']

    print(f2['velx'])
    velx = f2['velx']
    print(f2['vely'])
    vely = f2['vely']
    print(f2['velz'])
    velz = f2['velz']

    print() #whitespace
    #print(f2.attrs['Coordinates']).decode("utf-8")

    #for i in f2.attrs["VariableNames"]:
    #    print(i)

    #TODO: Find individual data keys and examine their data types.

    dict = {
        "coordinates" : coords,
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

print("\n Printing dictionary . . . . . . . . . . . . . . . . . . . .")

filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0000"
print(fileInputTest(filename))
filename = "/mnt/c/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/parkerCRs_hdf5_plt_cnt_0001"
print(fileInputTest(filename))
