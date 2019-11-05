#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get2dSliceNP.py
#  
#  Copyright 2019 Eric Eckert <dukeleto@spica>
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

"""
I have changed this version so that it gets the adjacent planes (in addition to just the main plane!)
"""

import sys
import h5py
import numpy as np
import datetime
import pandas as pd
import itertools 
import scipy
import os
from scipy import integrate

def dframeToDict(dFrame):
    #takes the input dataframe, unfucks it, and returns a dictionary
    dFrame = list(dFrame.iterrows())
    return {i[1].to_list()[0] : i[1].to_list()[1] for i in dFrame}
          
#stuff to insert into main

sw4essiout_filename = "Location14.cycle=0000.essi"
drm_filename = "DRM_slice_model_142D.hdf5"
parameterFile = "get2dSlice.csv"
#load up the parameter file
parameterFile = dframeToDict(pd.read_csv(parameterFile))
#check that the azimuth is valid
assert (int(parameterFile['azimuth'])==90 or int(parameterFile['azimuth'])==0)
if(int(parameterFile['azimuth'])==90):
    #I am using the sw4 x plane direction |
    directivity = 'sw4_x'
elif(int(parameterFile['azimuth'])==0):
    #I am using the sw4 y plane direction --
    directivity = 'sw4_y'

#load up the sw4 outputs
sw4essiout = h5py.File(sw4essiout_filename, 'r')
h  = sw4essiout['ESSI xyz grid spacing'][0]
x0 = sw4essiout['ESSI xyz origin'][0]
y0 = sw4essiout['ESSI xyz origin'][1]
z0 = sw4essiout['ESSI xyz origin'][2]
t0 = sw4essiout['time start'][0]
npts = sw4essiout['vel_0 ijk layout'].shape[0]
dt = sw4essiout['timestep'][0]
t1 = dt*(npts-1)
print ('grid spacing, h: ', h)
print ('ESSI origin x0, y0, z0: ', x0, y0, z0)
print('timing, t0, dt, npts, t1: ', t0, round(dt,6), npts, round(t1,6) )
print('Shape of HDF5 data: ', sw4essiout['vel_0 ijk layout'].shape)

nt = sw4essiout['vel_0 ijk layout'].shape[0]
time = np.linspace(t0, t1, npts)

#load the input data
plane = int(int(parameterFile['horizontal_plane'])/h)

if(directivity == 'sw4_x'):
    #n_horrizontal = sw4essiout['vel_0 ijk layout'].shape[1]
    #take all timesteps, and the plane as specified in the input file
    #get planes behind and infront of specified plane!
    data_horrizontal = sw4essiout['vel_0 ijk layout'][:,:,plane-2:plane+1,:]
    #set out of plane data to 0
    data_vertical = -1.0*sw4essiout['vel_2 ijk layout'][:,:,plane-2:plane+2,:]

else:
    #n_horrizontal = sw4essiout['vel_1 ijk layout'].shape[0]
    data_horrizontal = sw4essiout['vel_1 ijk layout'][:,plane-2:plane+1,:,:]
    #set out of plane data to 0
    data_vertical = -1.0*sw4essiout['vel_2 ijk layout'][:,plane-2:plane+1,:,:]

#check the spacing requested and the input spacing, if they disagree interpolate the inpute data!
if(int(parameterFile["target_h"]) != int(sw4essiout["ESSI xyz grid spacing"][0])):
    #interpolate with the specified algorithm
    print("Interpolating data from " + str(sw4essiout["ESSI xyz grid spacing"][0])
     + " to " + str(parameterFile["target_h"]) )
    
    #

#load up the geometery file
drm_file = h5py.File(drm_filename,"a")
coordinates = drm_file['Coordinates']
n_coord = int(coordinates.shape[0] / 3)
print('# of coordinates: ', n_coord)
#create data sets for outputs (as necessary)
if not "Time" in drm_file.keys():
    Time_dset          = drm_file.create_dataset("Time", data=time)
    Accelerations_dset = drm_file.create_dataset("Accelerations", (n_coord*3, nt))
    Displacements_dset = drm_file.create_dataset("Displacements", (n_coord*3, nt))
else:
    del drm_file['Time']
    Time_dset = drm_file.create_dataset("Time", data=time)
    del drm_file['Accelerations']
    Accelerations_dset = drm_file.create_dataset("Accelerations", (n_coord*3, nt))
    del drm_file['Displacements']
    Displacements_dset = drm_file.create_dataset("Displacements", (n_coord*3, nt))

#integrate all velocity data to form displacement data set
data_horrizontal_disp = scipy.integrate.cumtrapz(y=data_horrizontal[:,:,:,:],dx=dt,initial=0,axis=0)
data_vertical_disp = scipy.integrate.cumtrapz(y=data_vertical[:,:,:,:],dx=dt,initial=0,axis=0)
#differentiate all displacement data to form acceleration
data_horrizontal_acc = np.gradient(data_horrizontal[:,:,:,:],dt,axis=0)
data_vertical_acc = np.gradient(data_vertical[:,:,:,:],dt,axis=0)

"""
Check that the integration and differentiation look reasonalbe for a point
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.plot(time, data_horrizontal_disp[:,0,0,0])
fig.savefig("disp_plot_test.png")

fig, ax = plt.subplots()
ax.plot(time, data_horrizontal_acc[:,0,0,0])
fig.savefig("acc_plot_test.png")
"""
#save each coordinate to the container
coordinateList = coordinates[:]
coordinateList = list(coordinateList.reshape(coordinateList.shape[0]//3,-1,3))
#coordinateList = sorted(coordinateList,key=lambda coordinateList: coordinateList[0][2])   
for coordinate in range(len(coordinateList)): 
    #check which way is horrizontal
    if(directivity == 'sw4_x'):
        #fix the values if I have a negetive in my expression
        #checkValues = lambda x:0 if x != 0 else x         
        loc_i,loc_j,loc_k = int(coordinateList[coordinate][0][0]/h),int(coordinateList[coordinate][0][1]/h),int((int(parameterFile['essi_z_end'])-coordinateList[coordinate][0][2])/h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0 
        #note that j is unuesd !!
        #check my coordinate value  
        #loc j is out of plane
        plane = 1 + loc_j      
        Accelerations_dset[coordinate*3,:] =data_horrizontal_acc[:,loc_i,plane,loc_k] 
        Accelerations_dset[coordinate*3+1,:] = 0.0
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,loc_i,plane,loc_k]  
        
        Displacements_dset[coordinate*3,:] = data_horrizontal_disp[:,loc_i,plane,loc_k] 
        Displacements_dset[coordinate*3+1,:] = 0.0
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,loc_i,plane,loc_k] 
        print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_i,loc_k)) + "assinged to " + str(coordinate))
    
    else:
        #note that i is unuesd !!
        #checkValues = lambda x:0 if x != 0 else x         
        loc_j,loc_i,loc_k = int(coordinateList[coordinate][0][0]/h),int(coordinateList[coordinate][0][1]/h),int((int(parameterFile['essi_z_end'])-coordinateList[coordinate][0][2])/h)
        #z stop if z is ever negetive!
        #assert loc_k >= 0
        #get the correct plane based on the location coordinate
        #note that for this orientation loc_i is out of plane, for the other orientation
        #loc i is out of plane
        plane = 1 + loc_i
        Accelerations_dset[coordinate*3,:] = data_horrizontal_acc[:,plane,loc_j,loc_k]  
        Accelerations_dset[coordinate*3+1,:] = 0.0
        Accelerations_dset[coordinate*3+2,:] = data_vertical_acc[:,plane,loc_j,loc_k] 
        
        Displacements_dset[coordinate*3,:] = data_horrizontal_disp[:,plane,loc_j,loc_k]  
        Displacements_dset[coordinate*3+1,:] = 0.0
        Displacements_dset[coordinate*3+2,:] = data_vertical_disp[:,plane,loc_j,loc_k] 
        print("Coordinate pair " + str(coordinateList[coordinate][0]) + " computed at " + str((loc_j,loc_k)) + "assinged to " + str(coordinate))

        
#close the container
#sw4essiout.close()
#drm_file.close()


"""
def main(args):
    #load up the parameter csv file
    
    
    
    return 0

if __name__ == '__main__':
    #main()
    pass
"""