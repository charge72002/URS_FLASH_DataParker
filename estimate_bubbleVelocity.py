# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 10:53:46 2021

@author: wongb
"""

import math

def calcGravity(y):
    C = (-2 * math.pi) * (6.6743 * 10**-8) * (9.84 * 10**-3)
    H = 7.715 * 10**20
    return C * math.tanh(y/H)

def calcVelocity(y):
    g = calcGravity(y)
    v = math.sqrt( abs(2*y*g) )
    if(g<0): v*=-1;
    return v