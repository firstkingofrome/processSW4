#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  generate_boundary_text_files.py
#  
#  Copyright 2019 Eric Eckert 
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import re
import itertools
#import os

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

#fileName = "load.fei"
#section = "Boundary_surface"
def getsectionFEIfile(fileName,section):
    begunSection = False
    sectionData = []
    #load up the entire file to a list
    feiFile = open(fileName,'r')
    feiFile = feiFile.readlines()
    #get the current section that I am interested in
    for line in feiFile:
        if(not(begunSection) and section in line):
            begunSection = True
        elif(begunSection and section in line):
            begunSection = False
            break
        if(begunSection and "//" not in line):
            sectionData.append(line)
    sectionData = [i for i in sectionData if(len(i))>1]
    return sectionData


def listTotxtFile(fileName,lines):
    fileName = open(fileName,"w")
    fileName.write('\n'.join(lines))
    fileName.close()
    return

nodes = "node.fei"
node_coordinates = coordinatesToDict(nodes)
loadFEIFilename = "load.fei"
coordinates = []
# note that you can modify this if you would like to extract additional sections
#extractionSections specifies the section followed by its destination text tile while
#drm_elementSections specifies the DRM elements to extracted to the DRMelement.txt file!
#note that the order maters here! interior then exteriror
extractionSections = {"Boundary_surface":"DRMbound.txt","Exterior_surface":"DRMexterior.txt"}
drm_elementSections = {"DRMelement.txt":["#DRM1","#DRM2","#DRM3","#DRM4"]}
for keys in extractionSections:
    #get the data
    data = getsectionFEIfile(loadFEIFilename,keys)
    #get the coordiantes for this data
    data = [i.split(" ")[3] for i in data]
    #save the data to the specified text file
    listTotxtFile(extractionSections[keys],data)
    #get the node coordinates
    for node in data:
        coordinates.append(node_coordinates[node])
#now save the coordinates to the specified text file
coordinates = list(itertools.chain.from_iterable(coordinates))
#cast all as text
coordinates = [str(i) for i in coordinates]
listTotxtFile("coordinates.txt",coordinates)

elementData = []   
for keys in drm_elementSections:
    #loop through all of the elements specified in the drm elements list
    for element in drm_elementSections[keys]:
        elementData.append(getsectionFEIfile(loadFEIFilename,element))
    #flatten the lists to form one master list
    elementData = list(itertools.chain.from_iterable(elementData))
    #get only the integers
    elementData = [i.split(" ")[7].split(";")[0] for i in elementData]
    listTotxtFile(keys,elementData)

### now save all of the coordinates 
