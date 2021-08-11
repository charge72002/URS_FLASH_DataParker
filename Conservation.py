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
#%%s
## Mac filenames
filename = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0080"
filedirectory = "/Users/bwong/Downloads/URS_Data/m2_c1_16x8_64x64/More Plot Files"
pwd = "/Users/bwong/URS_FLASH_DataParker"
#%%
import yt
#adjust units as needed, otherwise it's alot
from yt.units import cm, kpc, erg, g, dyne
import os
import os.path
from os import path
import time
import matplotlib
import matplotlib.pyplot as plt
import beepy #sound for when the code is done running
import numpy as np
import math
import csv

import sys
sys.path.insert(0, pwd)
# Sherry packages
import estimate_bubbleVelocity as est
import hdf5_parser

#mess with classes
"""The energies class calculates simulation energy in ergs."""
class energies:
    #CONSTANTS
    #C= -2pi * newton's constant G * surface density \Sigma
    C = (-2 * math.pi) * (6.6743 * 10**-8) * (9.84 * 10**-3)
    H = 7.715 * 10**20
    
    ad, KE, PE_grav, PE_mag, PE_cr, PE_thermo = None, None, None, None, None, None
    """class constructor using yt dataset ds
    ~3-5 seconds to run per file"""
    def __init__(self, ds, gamma_cr = 4/3, gamma_gas = 5/3, C = C, H=H):
        startTime = time.time()
        # C = self.C
        # H = self.H
        # gamma_cr = 4/3
        # gamma_gas = 5/3
        ad = ds.all_data()
        
        #energy densities in dynes/cm^2 aka erg/cm^3
        KE = ad['kinetic_energy']
        magEnergy = ad[('gas', 'magnetic_energy')]
        pressure = ad[('gas', 'pressure')]
        #CR energy density in erg/g
        cray = ad['flash', 'cray'] * (erg/g)
        #gravitational potential in erg/g
        #use numpy instead of math package on array
        phi = C*H*np.log(np.cosh(ad['y']/H)) * erg/g
        
        volume = ad[('gas', 'cell_volume')]
        mass = ad[('gas', 'cell_mass')]
        
        
        #turn everything into ergs
        self.KE = KE * volume
        self.PE_grav = phi * mass
        # self.PE_mag =  ((magField**2) / (8*math.pi)) * volume #using Roark mag field in gauss
        self.PE_mag = magEnergy * volume
        # self.PE_cr = cray * (1/(gamma_cr-1)) * m #using Roark cray pres
        self.PE_cr = cray * mass
        self.PE_thermo = pressure * (1/(ds.gamma-1)) * volume
        print("energies.__init__() done. Time elapsed (sec): " + str(time.time()-startTime))

        
    """Returns a dictionary with field totals, including the total energy combined"""
    def totals(self):
        totals = {'total':0,
                  'KE': sum(self.KE), 
                  'PE_grav': sum(self.PE_grav), 
                  'PE_mag': sum(self.PE_mag), 
                  'PE_cr': sum(self.PE_cr), 
                  'PE_thermo': sum(self.PE_thermo)}
        for field in totals.keys(): 
            # print(field)
            totals['total'] += totals[field]
            # if not(field=='total'): totals['total'] += totals[field]
        return totals

# # directory = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/"
# # directory = "D://URS_LargeData/Parker_forSherry/"
directory = "D://URS_LargeData/m2_c1_16x8_64x64_cfl06_tx10/"
# os.path.exists(directory)
# # directory = "/Users/wongb/Documents/URS Data/diffusion_3e28/diffusion_3e28/"
# # directory = "/Users/wongb/Documents/URS Data/m1.5_c1_16x16_128x128_Rodrigues_Streaming/More Plot Files/"
saveDirectory = "D:/URS_LargeData/SherryPlots"
# path.exists(saveDirectory)

#%%
filenameTemp = "/Users/wongb/Documents/URS Data/m2_c1_16x8_64x64/More Plot Files/parkerCRs_hdf5_plt_cnt_0001"
ds = yt.load(filenameTemp)
ad = ds.all_data()

#%%
### test energies class
eng = energies(ds)
startTime = time.time()
totals = eng.totals()
print("energies.totals() done. Time elapsed (sec): " + str(time.time()-startTime))
# for field in totals.keys(): print(str(field) + ": " + str(totals[field]))
# print info on all fields
print("\n" + str(totals.keys()))
for field in totals.keys(): 
    print(str(field) + ": ")
    print("Total:\t" + str(totals[field]))
    if not(field == 'total'): 
        #vars() prints all variables in an object
        print("Max: \t" + str(np.max(vars(eng)[field])))
        print("Min: \t" + str(np.min(vars(eng)[field])))


#%%
### iterate total energies over time

startTime = time.time()
timeStamps = [] #list of time stamps
energylist = [] #list of energies objects
#dict of arrays of total energy over time
energytotals = {'total':[],
                'KE': [], 
                'PE_grav': [], 
                'PE_mag': [], 
                'PE_cr': [], 
                'PE_thermo': []}

#53 minute runtime. pain.
for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        timeStamps.append(int(timeStamp))
        ds = yt.load(filedirectory+ '/' +fileName)
        
        #that's a LOT of data (112 timesteps * 5 fields) please don't fry my computer
        #each file usually has 136 [derived] fields?
        #means ~43MB of data ok that's fine
        # energylist.append(energies(ds))
        
        temptotal = energies(ds).totals()
        for field in temptotal:
            energytotals[field].append(temptotal[field])
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Energy conservation done. Time elapsed (sec): " + str(time.time()-startTime))

#%%
### save important data as .csv
with open(pwd + '/Conservation/EnergiesTotal.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', "KE", "PE_grav", "PE_mag", "PE_cr", "PE_thermo", "total"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    #iterate through time steps
    for i in range(0, len(timeStamps)):
        #extract number w/o units
        def num(string): return float(str(string).split(" ")[0])
        eng = energytotals
        writer.writerow({'t':           timeStamps[i], 
                         "KE":          num(eng["KE"][i]),
                         "PE_grav":     num(eng["PE_grav"][i]),
                         "PE_mag":      num(eng["PE_mag"][i]),
                         "PE_cr":       num(eng["PE_cr"][i]),
                         "PE_thermo":   num(eng["PE_thermo"][i]),
                         "total":       num(eng["total"][i])    })
print(".csv output file written.")

#%%
### examine fields
for field in energytotals.keys(): 
    if(field != 't'):
        print(str(field) + ": ")
        #vars() prints all variables in an object
        print("\tAvg: \t" + str(np.mean(energytotals[field])))
        print("\tMax: \t" + str(np.max(energytotals[field])))
        print("\tMin: \t" + str(np.min(energytotals[field])))
        
#%%
## plot energy total
plt.close()
plt.plot(timeStamps, energytotals['total'])
plt.title("Total energy")
plt.savefig(pwd + '/Conservation/TotalEnergy')

#%%
## plot all energies on one plot
plt.close()
fig, ax = plt.subplots()
for field in energytotals:
    ax.plot(timeStamps, energytotals[field], label=field)
fig.suptitle("All energies")
ax.legend()
## fields in order
# ax.legend(['total', 'KE', 'PE_cr', 'PE_thermo','PE_mag', 'PE_grav'])
fig.savefig(pwd + '/Conservation/AllEnergies')
print("Plotting finished.")

#%%
## plot variation of energies from initial values
plt.clf()
plt.close()

initValues = {
            "KE":          energytotals["KE"][0],
            "PE_grav":     energytotals["PE_grav"][0],
            "PE_mag":      energytotals["PE_mag"][0],
            "PE_cr":       energytotals["PE_cr"][0],
            "PE_thermo":   energytotals["PE_thermo"][0],
            "total":       energytotals["total"][0]    
            }
variation = {'total':[],
            'KE': [], 
            'PE_grav': [], 
            'PE_mag': [], 
            'PE_cr': [], 
            'PE_thermo': []} 
for i in timeStamps:
    for field in variation.keys():
        variation[field].append( energytotals[field][int(i)] - initValues[field] )

fig, ax = plt.subplots()
for field in variation:
    ax.plot(timeStamps, variation[field], label=field)
# ax.set_yscale('symlog')
fig.suptitle("All energies; variation")
ax.legend()
fig.savefig(pwd + '/Conservation/VariationB')
print("Plotting finished.")


#%%
## just the total mass log for testing
fig, ax = plt.subplots()
ax.plot(timeStamps, variation['total'])
ax.set_yscale('symlog')
fig.suptitle("Total energy; log variation")

#%%
### READ csv
energytotals = {}
timestamps = []
with open(pwd + '/Conservation/EnergiesTotal.csv', 'r', newline='') as csvfile:
    r = csv.DictReader(csvfile)
    linenum=0
    for row in r:
        # print(row.keys())
        for label in row: 
            if(linenum==0): energytotals[label] = [] #init
            energytotals[label].append(float(row[label]))
        linenum+=1
    print(str(linenum) + " lines read.")
    print(energytotals.keys())
timeStamps = energytotals.pop('t')

#%%

################
# Examine fields
################

### compare derived KE field with classical KE calculation
# print(ds.field_list)
# print(ds.derived_field_list)

print(ad[('gas', 'kinetic_energy')][0] * ad[('gas', 'cell_volume')][0])
m = ad[('gas', 'cell_mass')][0]
# vol = ad[('gas', 'cell_volume')][0]
v = ad[('gas', 'velocity_magnitude')][0]
print(0.5 * m * (v**2))
# print(vol * 0.5 * m * (v**2))

#KE and classical KE nearly identical, therefore
#Thermal energy NOT INCLUDED in KE; must track separately
#%%
### repeat for full array

print("\nTOTAL energies:")
print(sum(ad[('gas', 'kinetic_energy')])* ad[('gas', 'cell_volume')][0])
m = ad[('gas', 'cell_mass')] #in g
v = ad[('gas', 'velocity_magnitude')] #in cm/s
print(sum(0.5 * m * (v**2))) #be sure to sum AFTER calculation

## access units
print("\nUnits: " + str(ds.mass_unit) + ", " + str(ds.length_unit) + ", " + str(ds.temperature_unit))
print('cray units in ergs/g')
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
# ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
# ad = ds.all_data()
for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(filedirectory+fileName)
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
# unit = str(ad[fields]).split("] ")[1]
startTime = time.time()
timeStamps = []
f = []
for field in fields: f.append([])
classic = []
ds = yt.load("D://URS_LargeData/Parker_forSherry/parkerCRs_hdf5_plt_cnt_0700")
ad = ds.all_data()
for fileName in os.listdir(filedirectory):
    if(fileName.startswith("parkerCRs")):
        #Start
        print(fileName)
        timeStamp = fileName[len(fileName)-4: len(fileName)]
        ds = yt.load(filedirectory+ '/' +fileName)
        ad = ds.all_data()
        timeStamps.append(int(timeStamp))
        for i in range(0, len(fields)):
            f[i].append(sum(ad[fields[i]]))
        #calculate classical energy
        m = ad[('gas', 'cell_mass')] #in g
        v = ad[('gas', 'velocity_magnitude')] #in cm/s
        classic.append(sum(0.5 * m * (v**2)))
beepy.beep(4)
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Energy conservation done. Time elapsed (sec): " + str(time.time()-startTime))

#%%
### save important data as .csv
with open(pwd + '/Conservation/Energies.csv', 'w', newline='') as csvfile:
    fieldnames = ['t', 'kinetic_energy', 'classical_energy', 'magnetic_energy']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for i in range(0, len(f[0])):
        writer.writerow({'t': i, 
                         'kinetic_energy': f[0][i],
                         'classical_energy': classic[i], 
                         'magnetic_energy': f[1][i]})

#%%
### READ csv
f = {}
with open(pwd + '/Conservation/Energies.csv', 'r', newline='') as csvfile:
    r = csv.DictReader(csvfile)
    linenum=0
    for row in r:
        for label in row: 
            if(linenum==0): 
                #init dict
                f[label] = []
            else:
                f[label].append(row[label])
        linenum+=1
    print(str(linenum) + " lines read.")
    print(f.keys())
    
#%%
### plot
plt.clf()
plt.plot(f['t'], f['kinetic_energy'], label = 'kinetic_energy', linestyle = '--', lw=3)
plt.plot(f['t'], f['classical_energy'], label = 'classical_energy', linestyle = ':', lw=3)
plt.plot(f['t'], f['magnetic_energy'], label = 'magnetic_energy')
plt.legend()

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

