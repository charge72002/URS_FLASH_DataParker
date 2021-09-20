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
sys.path.insert(0, pwd)

import athena_read as ath
import hdf5_parser


#%%

#
# data = ath.athdf(filename)
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
# ZorderX = out['posXarray']
# ZorderY = out['posYarray']
ZorderX = [0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3]
ZorderY = [0, 0, 1, 1, 0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 3, 3]

#setup large sequence
xlen = int( np.sqrt(len(ZorderX)) )
A = np.asarray(seqGenerator(xlen))
B = np.multiply(2, A)
twoD = []
for b in B:
    tempArray = [] 
    for a in A:
        tempArray.append(a+b)
    twoD.append(tempArray)
print(len(twoD))

#%%
coords = np.empty((xlen, xlen), dtype=float)
X = np.empty((xlen, xlen), dtype=float)
Y = np.empty((xlen, xlen), dtype=float)
for y in range(0, len(twoD)): #rows
    for x in range(0, len(twoD[y])): #cols
        index = twoD[y][x]
        xpos = ZorderX[index]
        ypos = ZorderY[index]
        
        X[y][x] = xpos
        Y[y][x] = ypos

    
#%%
plt.clf()
plt.plot(X, Y)
# thingX = [X[5][5],  X[6][5], X[5][6], X[6][6]]
# thingY = [Y[5][5],  Y[6][5], Y[5][6], Y[6][6]]
# thingX = [X[0][0],  X[0][1], X[1][0], X[1][1]]
# thingY = [Y[0][0],  Y[0][1], Y[1][0], Y[1][1]]
# plt.plot(thingX, thingY)

