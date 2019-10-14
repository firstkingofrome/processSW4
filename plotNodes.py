import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py

#nodes to plot
nodes = [[10,0,120],[10,-10,120],[10,10,120]]
#hdf5 file
geometeryFile = h5py.File("DRM_slice_model_142D.hdf5",'r')
coordinateList = geometeryFile["Coordinates"][:]
coordinateList = list(coordinateList.reshape(coordinateList.shape[0]//3,-1,3))
plotIndexes = []
for node in range(len(coordinateList )):
    #find the indexes to plot on
    if(list(coordinateList[node][0]) in nodes):
        plotIndexes.append(node)
        

#plot accelerations
fig, ax = plt.subplots()
ax=plt.gca()
ax.plot(geometeryFile["Time"][:],geometeryFile['Accelerations'][plotIndexes[0]*3],label=str(nodes[0]))
ax.plot(geometeryFile["Time"][:],geometeryFile['Accelerations'][plotIndexes[1]*3],label=str(nodes[1]))
ax.plot(geometeryFile["Time"][:],geometeryFile['Accelerations'][plotIndexes[2]*3],label=str(nodes[2]))
fig.show()
fig.savefig("acc_plot_test.png")


#plot displacement
fig, ax = plt.subplots()
ax=plt.gca()
ax.plot(geometeryFile["Time"][:],geometeryFile['Displacements'][plotIndexes[0]*3],label=str(nodes[0]))
ax.plot(geometeryFile["Time"][:],geometeryFile['Displacements'][plotIndexes[1]*3],label=str(nodes[1]))
ax.plot(geometeryFile["Time"][:],geometeryFile['Displacements'][plotIndexes[2]*3],label=str(nodes[2]))
fig.show()
fig.savefig("acc_plot_test.png")


geometeryFile.close()
#plot displacements


#now plot all of the stuff in those indexes

