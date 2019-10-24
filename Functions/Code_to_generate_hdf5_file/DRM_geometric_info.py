# -*- coding: utf-8 -*-
# @Author: Suiwen Wu
# @Date:   2019-04-28 11:46 pm


#! /usr/bin/env python

import h5py
import scipy as sp
import time
#import pickle
#import sys
#import math

#from obspy import read


## Write elements and ndoes data
elements = sp.loadtxt("DRMelements.txt",dtype=sp.int32)
exterior_nodes = sp.loadtxt("DRMexterior.txt",dtype=sp.int32)
boundary_nodes = sp.loadtxt("DRMbound.txt",dtype=sp.int32)
Coordinates = sp.loadtxt("Coordinates.txt",dtype=sp.float_)

Ne = sp.array(exterior_nodes.size)
Nb = sp.array(boundary_nodes.size)

Nt =Ne+Nb

all_nodes = sp.hstack((boundary_nodes, exterior_nodes))
is_boundary_node = sp.zeros(Nt,dtype=sp.int32)
is_boundary_node[0:Nb]=1

h5file=h5py.File("DRM_slice_model_14.hdf5","w")

h5file.create_dataset("Coordinates",data=Coordinates)
h5file.create_dataset("DRM Nodes",data=all_nodes)
h5file.create_dataset("Elements",data=elements)


h5file.create_dataset("Is Boundary Node",data=is_boundary_node)

h5file.create_dataset("Number of Exterior Nodes",data=Ne)
h5file.create_dataset("Number of Boundary Nodes",data=Nb)