import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import datetime
import pandas as pd
import scipy
from scipy import integrate


def interpolate(val0, val1, fraction):
    return val0 + (val1 - val0) * fraction

sw4essiout_filename = "Location15.cycle=00000.essi"
sw4essiout = h5py.File(sw4essiout_filename, 'r')
t0 = sw4essiout['time start'][0]
npts = sw4essiout['vel_0 ijk layout'].shape[0]
dt = sw4essiout['timestep'][0]
t1 = dt*(npts-1)
time = np.linspace(t0, t1, npts)
nt = sw4essiout['vel_0 ijk layout'].shape[3]
dataVertical = -1.0*sw4essiout['vel_2 ijk layout'][:,:,10,:]
dataHorizontal =sw4essiout['vel_0 ijk layout'][:,:,10,:]

###compute with numpy
dispVerticalNumpy = np.copy(dataVertical)
dispVerticalNumpy = scipy.integrate.cumtrapz(y=dispVerticalNumpy[:,:,:],dx=dt,initial=0,axis=0)

accVerticalNumpy = np.copy(dataVertical)
accVerticalNumpy = np.gradient(accVerticalNumpy[:,:,:],dt,axis=0)
#accVerticalNumpy = np.diff(accVerticalNumpy[:,:,:],n=1,axis=1)

###compute manually
dispVerticalNoNumpy = np.copy(dataVertical) 
accVerticalNoNumpy = np.copy(dataVertical) 
### manually compute difference and integrand
prev_disp = 0
current_disp = 0
prev_accel = 0
current_accel = 0
for t in range(1,nt-1):
    plane_vel = dispVerticalNoNumpy[t,0,0]
    prev_plane_vel = dispVerticalNoNumpy[t-1,0,0]
    next_plane_vel = dispVerticalNoNumpy[t+1,0,0]
    if t > 1:
        prev_disp = current_disp   
    current_disp = prev_disp + 0.5 * (prev_plane_vel + next_plane_vel) * dt
    dispVerticalNoNumpy[t,0,0] = current_disp
    current_accel = (next_plane_vel - prev_plane_vel) / (2.0 * dt)
    accVerticalNoNumpy[t,0,0]= current_accel


fig, ax = plt.subplots()
ax.plot(time, dispVerticalNoNumpy[:,0,0],label="numpy data")
ax.plot(time, dispVerticalNumpy[:,0,0],label="non numpy data (1st order)")
ax.legend(loc='upper left')
fig.savefig("disp_comparing_numpy.png")


fig, ax = plt.subplots()
ax.plot(time[0:10633], accVerticalNumpy[:,0,0][0:10633],color='r',label="numpy data")
ax.plot(time[0:10633], accVerticalNoNumpy[:,0,0][0:10633],color='b',label="non numpy data (1st order)")
ax.legend(loc='upper left')
fig.savefig("acc_comparing_numpy.png")


"""



"""

