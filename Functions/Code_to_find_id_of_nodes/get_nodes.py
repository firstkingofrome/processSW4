import sys
import h5py
import numpy as np
import datetime
import pandas as pd
import itertools 
import scipy
import re
from scipy import integrate
def dframeToDict(dFrame):
    #takes the input dataframe, unfucks it, and returns a dictionary
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list() for i in dFrame}

def coordinatesToDict(nodeFileName):        
    coordinates = open(nodeFileName,'r')
    coordinates = coordinates.readlines()
    coordinates = [i for i in coordinates if "//" not in i]
    #ditch every empty line
    coordinates = [i for i in coordinates if len(i)>1]
    #refromat as a dictionary of lists
    coordinateDict = {}
    for i in coordinates:
        nodeNumber = i.split("#")[1].split(" ")[1]
        coordinateDict[nodeNumber] = re.findall('([-^\d*]\d*[*])',i)
        coordinateDict[nodeNumber] = [int(index.strip("*")) for index in coordinateDict[nodeNumber]]
    return coordinateDict

#gets the node number based on the xyz coordinates and returns it
def getNodeNumber(nodes,x,y,z):
    for key in nodes.keys():
        if(nodes[key][0]==x and nodes[key][1]==y and nodes[key][2]==z):
            print(key)
            return int(key)
    print("WARNING the node coordinate specified was missing from the node.fei master file")
    print(x,y,z)
    return "NoneInMasterFile"
        
parameterFileName = "Comparison_points.csv"
nodes = "node.fei"

### load the parameter file
parameterFile = dframeToDict(pd.read_csv(parameterFileName))
### Load the .feioutput filecoordinates 
nodes = coordinatesToDict(nodes)
### asscoiate a coordinate with a node number
nodeNumbers = {}
for line in parameterFile.keys():
    x,y,z = parameterFile[line][1:4]
    nodeNumbers[line] = getNodeNumber(nodes,x,y,z)

#now write all of the specified nodes in the master file
#open the master csv file for writing, this is very inelligent, but what can I say?
csvFileEditable = pd.read_csv(parameterFileName)
for index,row in csvFileEditable.iterrows():
    csvFileEditable.loc[index,"Node_#"] = nodeNumbers[csvFileEditable.loc[index,"point"]]
    
#now save the updated csv file
csvFileEditable.to_csv(parameterFileName,index=False)
