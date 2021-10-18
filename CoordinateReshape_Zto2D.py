#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 15:12:28 2021

@author: bwong
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

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker, cm #for log scale
sys.path.insert(0, pwd)
# import copy

import athena_read as ath
import h5py
import hdf5_parser


#%%

#
# data = ath.athdf(filename)
rawhdf = h5py.File(filename, 'r')
out = hdf5_parser.setup(filename)

#%%
### failed reshaping attempts
X = np.reshape(out['posXarray'], [512, 512], order ='C')
Y = np.reshape(out['posYarray'], [512, 512], order ='C')

thingX = [X[5, 5],  X[6, 5], X[5, 6], X[4, 5], X[5, 4]]
thingY = [Y[5, 5],  Y[6, 5], Y[5, 6], Y[4, 5], Y[5, 4]]


plt.plot(X, Y)
# plt.plot(thingX, thingY)

#%%
### Math the reshaping
### Moser-De Brujin sequence
### see https://www.wikiwand.com/en/Moser–De_Bruijn_sequence

#Generates a Moser-De Brujin sequence w/ length n
def seqGenerator (n):
    result = [0, 1]
    for i in range(2, n):
        if i%2==0:  #even
            # print("even " + str(int(i/2)))
            result.append(4*result[int(i/2)])
        else:       #odd
            # print("odd " + str(int((i/2))-1))
            result.append(4*result[int((i/2))]+1)
    return result

A = np.asarray(seqGenerator(8))
B = np.multiply(2, A)
twoD = []
for b in B:
    tempArray = [] 
    for a in A:
        tempArray.append(a+b)
    twoD.append(tempArray)
#IT WORKS MWAHAHAHA

#%%
### New Z order to cartesian attempt

#test with small array first
# ZorderX = [0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3]
# ZorderY = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3]
# parsed array from hdf5_parser
# ZorderX = out['posXarray']
# ZorderY = out['posYarray']
# raw hdf5 array
ZorderX = []
ZorderY = []
for xyz in rawhdf['coordinates']:
    ZorderX.append(xyz[0])
    ZorderY.append(xyz[1])


#setup large sequence
xlen = int( np.sqrt(len(ZorderX)) ) #ASSUME SQUARE BOX
A = np.asarray(seqGenerator(xlen))
B = np.multiply(2, A)
twoD = []
for b in B:
    tempArray = [] 
    for a in A:
        tempArray.append(a+b)
    twoD.append(tempArray)
print(len(twoD))

#assign Z indicies to coordinate array
coords = np.empty((xlen, xlen), dtype=float)
origX = np.empty((xlen, xlen), dtype=float)
origY = np.empty((xlen, xlen), dtype=float)
for y in range(0, len(twoD)): #rows
    for x in range(0, len(twoD[y])): #cols
        index = twoD[y][x]
        xpos = ZorderX[index]
        ypos = ZorderY[index]
        origX[y][x] = xpos
        origY[y][x] = ypos


## Create 8x8 full resolution afterwards
stepX = origX[0][1] - origX[0][0]
stepY = origY[1][0] - origY[0][0]
# print("X step: " + str(stepX))
# print("Y step: " + str(stepY)) #step before 8x8 splitup
X = np.empty((xlen*8, xlen*8), dtype=float)
Y = np.empty((xlen*8, xlen*8), dtype=float)


arrayin = origX
ylen = len(arrayin)
xlen = len(arrayin[0])
# out = np.empty((ylen, xlen), dtype=float)

for row in range(0, len(origX)):
    for col in range(0, len(origX[row])):
        tempXlin = np.linspace(origX[row][col], origX[row][col] + stepX, 8)
        tempYlin = np.linspace(origY[row][col], origY[row][col] + stepY, 8)
        tempMeshgrid = np.meshgrid(tempXlin, tempYlin)
        #iterate row/col of small meshgrid
        for subrow in range(0, len(tempMeshgrid[0])):
            for subcol in range(0, len(tempMeshgrid[0][0])):
                #Transform origX coords into X coords (8x8 bigger)
                ROW = (row*8)+subrow
                COL = (col*8)+subcol
                X[ROW][COL] = tempMeshgrid[0][subrow][subcol]
                Y[ROW][COL] = tempMeshgrid[1][subrow][subcol]
X, Y = (np.asarray(X), np.asarray(Y))

    

plt.clf()
plt.plot(X, Y)
# thingX = [X[5][5],  X[6][5], X[5][6], X[6][6]]
# thingY = [Y[5][5],  Y[6][5], Y[5][6], Y[6][6]]
# thingX = [X[0][0],  X[0][1], X[1][0], X[1][1]]
# thingY = [Y[0][0],  Y[0][1], Y[1][0], Y[1][1]]
# plt.plot(thingX, thingY)

##RESULT? No good. Z order combined with existing hdf5_parser code
##        mixed up 8x8 cartesian and z order arrays.
#%%
### test defined method

X, Y = hdf5_parser.cartesianCoordinates(h5py.File(filename, 'r'))
plt.clf()
# plt.plot(X, Y)
# colormap = np.linspace(1, len(X)*512, len(X)*512)
#repeat 2D 512 times
colormap = np.meshgrid(np.linspace(1, len(X), len(X)), np.empty(len(X), dtype=float))[0]
plt.scatter(X, Y, s=1, c=colormap)
plt.savefig(pwd + "/Plots/Zclean2.png")

#%%
newOut = hdf5_parser.setup(filename, format = "cartesian")
x = newOut['posXarray'].flatten()
y = newOut['posYarray'].flatten()
z = newOut['densityArray'].flatten()
# from matplotlib import ticker, cm #for log scale
# plt.tricontourf(x, y, z, levels=100, locator=ticker.LogLocator())

lev = np.logspace(np.log10(z.min()), np.log10(z.max()), num=1000)
plt.clf()
plt.tricontourf(x, y, z, locator=ticker.LogLocator(), levels = lev) #good for irregular Z values
plt.savefig(pwd + "/Plots/ZcleanB1.png")
