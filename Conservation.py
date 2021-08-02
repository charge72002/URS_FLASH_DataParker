# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 14:10:46 2021

@author: SherryWong
"""
#%%
## Windows Filenames
filename = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
#no slash at end
filedirectory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/wongb/Documents/Python_Scripts"
#%%
## Mac filenames
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/bwong/URS_FLASH_DataParker"
#%%
import yt
import os
import os.path
from os import path
import time
import matplotlib
import matplotlib.pyplot as plt
import beepy #sound for when the code is done running
import numpy as np
import math

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import hdf5_parser

# directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
# directory = "D://URS_LargeData/Parker_forSherry/"
directory = "D://URS_LargeData/m2_c1_16x8_64x64_cfl06_tx10/"
os.path.exists(directory)
# directory = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/"
# directory = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
path.exists(saveDirectory)

#%%
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()

#%%
# print(ds.field_list)
# print(ds.derived_field_list)

print(ad[('gas', 'kinetic_energy')][0])
m = ad[('gas', 'cell_mass')][0]
v = ad[('gas', 'velocity_magnitude')][0]
print(0.5 * m * (v**2))

print("\nTOTAL energies:")
print(sum(ad[('gas', 'kinetic_energy')]))
# m = sum(ad[('gas', 'cell_mass')]) #in g
# v = sum(ad[('gas', 'velocity_magnitude')]) #in cm/s
# print(0.5 * m * (v**2))
m = ad[('gas', 'cell_mass')] #in g
v = ad[('gas', 'velocity_magnitude')] #in cm/s
print(sum(0.5 * m * (v**2)))
#%%

#####################
#Quantity f over time
#####################
#this takes ~12 mins to run
field = ('gas', 'kinetic_energy')
unit = str(ad[field]).split("] ")[1]
startTime = time.time()
timeStamps = []
f = []
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        timeStamps.append(int(timeStamp))
        f.append(sum(ad[field]))
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Mass calc done. Time elapsed (sec): " + str(time.time()-startTime))

plt.clf()
# upperTick = round( np.amax(mass) )
# bottomTick = round( np.amin(mass) )
# log = np.logspace(bottomTick, upperTick, 8)
# plt.plot(timeStamps, np.log10(mass))
plt.plot(timeStamps, f)
# plt.semilogy(timeStamps, mass)
# plt.set_yscale('log', base=10)

plt.title("('gas', 'kinetic_energy')")
plt.xlabel("Timestamp")
# plt.ylabel("Total mass (g)")
plt.ylabel("Kinetic Energy (dynes/cm^2)")
# plt.ticklabel_format(axis="y", style="sci", scilimits=(38, 42), useMathText=True)
# plt.ticklabel_format(axis="x", style="plain")
# plt.yticks()
# fileName = "dm1.5_c1_16x16_128x128_Rodrigues_Streaming"
fileName = "m2_c1_16x8_64x64"
plt.savefig(pwd + "/Conservation/KE_"+fileName+".png")
plt.show()

#%%
## get everything. sum all energy-related terms over time, so you can manipulate them later

fields = [('gas', 'kinetic_energy'), ('gas', 'magnetic_energy')]
unit = str(ad[field]).split("] ")[1]
startTime = time.time()
timeStamps = []
f = np.array(len(fields))
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()
for fileName in os.listdir(directory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(directory+fileName)
        ad = ds.all_data()
        timeStamps.append(int(timeStamp))
        for i in range(0, len(fields)-1):
            f[i].append(sum(ad[fields[i]]))
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Mass calc done. Time elapsed (sec): " + str(time.time()-startTime))

#%%
### save important data as .csv
with open(pwd + '/Conservation/Energies.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', 'kinetic_energy', 'Ypos', 'prominences']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(0, len(total_peaks)):
        writer.writerow({'t': i+65, 
                         '#peaks': total_peaks[i],
                         'Ypos': Ypos[i], 
                         'prominences': prominences[i]})

#%%
###############################
#Calculate total mass SELECTION
###############################

startTime = time.time()
timeStamps = []
mass = []

out = hdf5_parser.setup(filename)
ymax = max(out['posYarray'])
ymin = min(out['posYarray'])
TFselect = np.logical_and(out['posYarray'] != ymax, out['posYarray'] != ymin)

for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        timeStamps.append(int(timeStamp))
        out = hdf5_parser.setup(filedirectory + "/" + fileName)
        ds = yt.load(filedirectory + "/" + fileName)
        ad = ds.all_data()
        mass.append(sum(ad[('gas', 'cell_mass')][TFselect]))
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Mass calc done. Time elapsed (sec): " + str(time.time()-startTime))
beepy.beep(4)

plt.clf()
plt.semilogy(timeStamps, mass)

plt.xlabel("Timestamp")
plt.ylabel("Total mass (g)")
fileName = "m2_c1_16x8_64x64"
plt.savefig(saveDirectory + "/totalMassSelect_m2_c1_16x8_64x64.png")
plt.show()

