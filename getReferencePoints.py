#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:48:57 2019
@author: eeckert
"""
### not flipped

import h5py
import numpy as np
import re
import pandas

sw4essiout_filename =  "Location5.cycle=00000.essi"
pointsRequested = "outputLocations.csv"
pointsOutput = "outputPoints"
#load the csv file
points = pandas.read_csv(pointsRequested)
#open the hdf5 container for this run and get each out put velocity point

sw4essiout = h5py.File(sw4essiout_filename, 'r')
spacing = int(sw4essiout['ESSI xyz grid spacing'][:][0])
#get constant values
sw4Xstart,sw4Ystart,sw4Zstart = int(points["sw4Xstart"][0]),int(points["sw4Ystart"][0]),int(points["sw4Zstart"][0])
openSeesxStart,openSeesyStart,openSeeszStart = int(points["openSeesxStart"][0]),int(points["openSeesyStart"][0]),int(points["openSeeszstart"][0])
for index, row in points.iterrows():
    print(row)
    #translate into xyz coordinates
    ### THE ESSI CORDINATES ARE XYZ WHILE SW4 IS YXZ
    OSx,OSy,OSz = row["x"],row["y"],row["z"]
    print("preparing refrence point " + row["point"])
    print("saving x y and z velocities for  exssi point " + "%i,%i,%i" % (OSx,OSy,OSz))
    """
    THE LINE BELOW THIS IS WHERE I CONVERT FROM ESSI --> SW4
    """
    sw4x,sw4y,sw4z = (sw4Xstart+(OSy-openSeesyStart))/spacing,(sw4Ystart+(OSx-openSeesxStart))/spacing,(sw4Zstart+(openSeeszStart-OSz))/spacing
    print("This corresponds to sw4 index %i,%i,%i"  %(sw4x,sw4y,sw4z))
    print('\n')
    #save xyz as csv files for that point
    #save ESSI x y and z!
    fname = "_vel_sw4_"+row["point"]+".csv"
    x = sw4essiout['vel_1 ijk layout'][:,sw4x,sw4y,sw4z]
    np.savetxt(fname='x'+fname,X=x)
    #save y
    y = sw4essiout['vel_0 ijk layout'][:,sw4x,sw4y,sw4z]
    np.savetxt(fname='y'+fname,X=y)
    #save z
    z = -1.0*sw4essiout['vel_2 ijk layout'][:,sw4x,sw4y,sw4z]
    np.savetxt(fname='z'+fname,X=z)
