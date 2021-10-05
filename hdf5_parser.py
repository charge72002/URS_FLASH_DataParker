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
    
def setup(fileName, format="Z"):
    file = h5py.File(fileName, 'r')
        
    print("Keys:")
    print(file.keys())
        
    if(format=="Z"):
        posXarray, posYarray = bigCoordinateSetup(file)
        
        return {"posXarray" : posXarray, 
                "posYarray" : posYarray, 
                "velXarray" : fieldParse(file, 'velx'), 
                "velYarray" : fieldParse(file, 'vely'), 
                "densityArray" : fieldParse(file, 'dens'), 
                "tempArray" : fieldParse(file, 'temp'),
                "magXarray" : fieldParse(file, 'magx'),
                "magYarray" : fieldParse(file, 'magy')}
    
    elif(format=="cartesian"):
        posXarray, posYarray = cartesianCoordinates(file)
        
        return {"posXarray" : posXarray, 
                "posYarray" : posYarray, 
                "velXarray" : ZtoCartesian(fieldParse(file, 'velx')), 
                "velYarray" : ZtoCartesian(fieldParse(file, 'vely')), 
                "densityArray" : ZtoCartesian(fieldParse(file, 'dens')), 
                "tempArray" : ZtoCartesian(fieldParse(file, 'temp')),
                "magXarray" : ZtoCartesian(fieldParse(file, 'magx')),
                "magYarray" : ZtoCartesian(fieldParse(file, 'magy'))}
    else:
        raise Exception("Invalid format " + format)
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

## For below methods, @see CoordinateReshapeAto2D.py
#helper method for ZtoCartesian() and cartesianCoordinates()
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

#Use Moser-De Brujin sequence to transform a 1D Z ordered array into a 
#2D cartesian array.
#takes a 1D Z-ordered array
#returns a 2D cartesian array
def ZtoCartesian(field):    
    #setup large sequence
    xlen = int( np.sqrt(len(field)) ) #ASSUME SQUARE BOX
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
    result = np.empty((xlen, xlen), dtype=float)
    for y in range(0, len(twoD)): #rows
        for x in range(0, len(twoD[y])): #cols
            index = twoD[y][x]
            result[y][x] = field[index]
    return np.asarray(result)
    

    return result
    
#apply method to coordinates in the file.
#accepts an hdf5 file
#returns 2D arrays X, Y
def cartesianCoordinates(file):
    ZorderX = []
    ZorderY = []
    #iterate and separate XYZ coordinates
    for xyz in file['coordinates']:
        ZorderX.append(xyz[0])
        ZorderY.append(xyz[1])
    xlen = int( np.sqrt(len(ZorderX)) ) #ASSUME SQUARE BOX
    
    origX = ZtoCartesian(ZorderX)
    origY = ZtoCartesian(ZorderY)
    
    # Create 8x8 full resolution afterwards
    stepX = origX[0][1] - origX[0][0]
    stepY = origY[1][0] - origY[0][0]
    X = np.empty((xlen*8, xlen*8), dtype=float)
    Y = np.empty((xlen*8, xlen*8), dtype=float)
    
    arrayin = origX
    ylen = len(arrayin)
    xlen = len(arrayin[0])
    
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
    return X, Y


