#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 10:08:47 2019
Clips the data set to the specified time range
@author: eeckert
"""
import sys
import h5py
import numpy as np
import datetime
import math


TIME_START=3.0
TIME_STOP=15.0

drm_filename = "Model_5m_mesh.hdf5"
drm_file = h5py.File(drm_filename,'r+')
#load the data sets
Time_dset          = drm_file['Time'][:]
Accelerations_dset = drm_file['Accelerations'][:]
Displacements_dset = drm_file['Displacements'][:]
#trim off the time
dt = Time_dset[-1]/(Time_dset.shape[0])
startStep = int(TIME_START/dt)
stopStep = int(TIME_STOP/dt)
#recompute time range
#Time_dset = Time_dset[startStep:stopStep] 

Accelerations_dset = Accelerations_dset[:,startStep:stopStep] 
Displacements_dset = Displacements_dset[:,startStep:stopStep] 
#save it!
del drm_file['Time']
del drm_file['Accelerations']
del drm_file['Displacements']


drm_file.create_dataset("Accelerations", data=Accelerations_dset)
drm_file.create_dataset("Displacements", data=Displacements_dset)

#generate the new time range
t0=0
t1=drm_file['Accelerations'].shape[1]*dt
time = np.linspace(t0,t1,drm_file['Accelerations'].shape[1])
#assign the new time data set
drm_file.create_dataset("Time",data=time)
print("Time Step="+str(dt))
print("Num steps="+str(Accelerations_dset.shape[0]))
drm_file.close()