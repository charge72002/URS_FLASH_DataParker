# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 00:03:42 2021

@author: wongb
"""

import numpy as np

"""check validity of LHS of gravitational potential equation
given a list of y positions and velocities, return an array of """
def LHS (y, v=0): #array-likes
    C = (-2 * np.pi) * (6.6743 * 10**-8) * (9.84 * 10**-3)
    H = 7.715 * 10**20
    y = np.asarray(y)
    
    #no velocity array passed; calculate here
    if(v==0):
        templist = list()
        convert = 3.1709792E-15 #cm/10Myr --> cm/s
        for i in range(1, len(y)):
            templist.append( (y[i-1] - y[i])*convert )
        v = np.asarray(templist)
        print(v)
    
    newy = list()
    vsquared = list()
    if(len(y) == (len(v)+1)): 
        for i in range(0, len(y)-1):
            newy.append((y[i] + y[i+1])/2)
            vsquared.append( v[i]**2 )
        # print(newy)
        # print(vsquared)
        newy = np.asarray(newy)
        return 0.5*(vsquared - (2*C*H*np.log(np.cosh(newy/H))))
    #lambda expression! 
    elif (len(y)==len(v)): vsquared = np.array(list(map(lambda x:pow(x, 2), v)))
    elif(len(y)!=len(v)): raise Exception("Incompatible array-likes with lengths (" + str(len(y)) + ", " + str(len(v)) + ")")
    # print(vsquared)
    return 0.5*(vsquared - (2*C*H*np.log(np.cosh(y/H))))

#upper bubble
y = [
-6.52808E+21,
-5.59156E+21,
-4.98558E+21,
-4.04906E+21,
-3.22272E+21]
print("\nUpper bubble:\n" + str(y))
print(LHS(y))

#lower bubble 1 (red)
y = [
-9.88853E+21,
-1.01640E+22,
-1.06047E+22,
-1.10454E+22]
print("\nLower bubble 1 (red):\n" + str(y))
print(LHS(y))

#lower bubble 2 (green)
y = [
-8.84183E+21,
-9.32520E+21,
-5.59156E+21]
print("\nLower bubble 2 (green):\n" + str(y))
print(LHS(y))