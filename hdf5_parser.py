# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 10:58:57 2021

@author: wongb
"""
import h5py
import numpy as np   
    
# coordinates get their own special setup to correct for their wonky shape
def bigCoordinateSetup(file):
    coordField = file['coordinates']
    #Make lists of individual data columns
    posX = [] #list
    posY = [] #list
    
    #xOffset = int()
    stepX = coordField[1][0] - coordField[0][0]
    stepY = coordField[2][1] - coordField[1][1] #bug was here
    print("X step: " + str(stepX))
    print("Y step: " + str(stepY))
    
    for coord in coordField:
        tempXlin = np.linspace(coord[0], coord[0] + stepX, 8)
        tempYlin = np.linspace(coord[1], coord[1] + stepY, 8)
        tempMeshgrid = np.meshgrid(tempXlin, tempYlin)
        posX.append(tempMeshgrid[0])
        posY.append(tempMeshgrid[1]) #bug is not here
        #print(coord)
            
    # print("Length of posX: " + str(len(posX)))
    posXarray, posYarray = (np.asarray(posX), np.asarray(posY))
    # print("Length of posX flattened: " + str(len(posXarray.flatten())))
    posXarray = posXarray.flatten()
    posYarray = posYarray.flatten()
    return posXarray, posYarray
    # END OF METHOD
    

# file is an hdf5 file.
# fieldname is some string that has a key in the data set
def fieldParse(file, fieldname):
    if not fieldname in file.keys():
        raise Exception("Key \"" + fieldname + "\" not found!" + "\nAcceptable keys:\n" + file.keys())
        
    fieldList = []
    for item in file[fieldname]:
        fieldList.append(item)
    return np.asarray(fieldList).flatten();
    
def setup(fileName):
    file = h5py.File(fileName, 'r')

    print("Keys:")
    print(file.keys())
    
    posXarray, posYarray = bigCoordinateSetup(file)
    
    return {"posXarray" : posXarray, 
            "posYarray" : posYarray, 
            "velXarray" : fieldParse(file, 'velx'), 
            "velYarray" : fieldParse(file, 'vely'), 
            "densityArray" : fieldParse(file, 'dens'), 
            "tempArray" : fieldParse(file, 'temp'),
            "magXarray" : fieldParse(file, 'magx'),
            "magYarray" : fieldParse(file, 'magy')}
    #END OF METHOD
    
# for a custom selection of fields
# comes with XY coordinates
def customSetup(fileName, fieldList):
    file = h5py.File(fileName, 'r')

    print("Keys:")
    print(file.keys())
    
    posXarray, posYarray = bigCoordinateSetup()
    
    result = {"posXarray" : posXarray, 
            "posYarray" : posYarray}
    
    for fieldName in fieldList:
        result[fieldName] = fieldParse(fieldName)
    
    return result;
    #END OF METHOD