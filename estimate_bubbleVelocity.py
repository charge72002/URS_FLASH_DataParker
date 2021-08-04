# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 10:53:46 2021

@author: wongb
"""
import numpy as np
import math

#CONSTANTS
#C= -2pi * newton's constant G * surface density \Sigma
C = (-2 * math.pi) * (6.6743 * 10**-8) * (9.84 * 10**-3)
H = 7.715 * 10**20

# returns gravitational acceleraton in cm/s^2
def calcGravity(y):
    return C * math.tanh(y/H)

# TODO: implement additional gravity within 1kpc (3E21cm) of galactic disk
# returns estimated ballistic velocity in cm/s
def calcVelocity(y, y0=0, v0=0):
    #gravity becomes constant after 4 scale heights
    if(abs(y)<(4*H)): #y<1~kpc for Solar Neighborhood Parameters
        print("Inside the disk 2")
        b = 2*C*H*math.log(math.cosh(y0/H))
        c = 2*C*H*math.log(math.cosh(y/H))
        X = v0**2 - b + c
        # print(X)
        #resolve sign issue
        if(X>0): return math.sqrt(X)
        else: return -(math.sqrt(-X))
    else:
        # g = calcGravity(y) 
        g = C #Constant gravity; tanh(y/H)~1
        v = math.sqrt( abs(2*y*g) )
        if(g<0): v*=-1;
        return v

# useful misc methods from BubbleVelocityInteractive.py
def removeDuplicates(arrayIN):
    listOUT = []
    for i in arrayIN:
        if(listOUT.count(i)==0):
            listOUT.append(i)        
    arrayOUT = np.asarray(listOUT)
    return arrayOUT
    #END OF METHOD

#Finds the number in the array closest to the target number.
def closestNum(arrayIN, targetNum):
    smallArray = removeDuplicates(arrayIN)
    currentDev = float("inf")
    result = float("inf")
    #TODO: if you wanna have fun make this a binary search thing
    for i in smallArray:
        newDev = abs(targetNum-i)
        if(newDev < currentDev):
            currentDev =  newDev
            result = i
    if(result == float("inf")):
        raise Exception("No match found")
    return result
    #END OF METHOD

## returns the indices right before a sign flip
## i.e. findSignFlips([-3, -1, 1, 3, -2]) returns [1, 3]
def findSignFlips(arrayIN):
    result = []
    prevNum = arrayIN[0]
    for index in range(0, len(arrayIN)-1):
       if prevNum < 0 and arrayIN[index] > 0: result.append(index-1)
       elif prevNum > 0 and arrayIN[index] < 0: result.append(index-1)
       prevNum = arrayIN[index]
    return np.array(result, dtype=int)