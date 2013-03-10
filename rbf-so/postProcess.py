# -*- coding: utf-8 -*-
import cPickle

from numpy import *
from random import *
from scipy import sparse,linalg
import os, sys
import matplotlib.pyplot as plt

#Enter which directory to post process
nPoints = 2450



fname = open('data_co2.pkl', 'r')
data = cPickle.load(fname)
fname.close()

fname = open('./Data/'+str(nPoints)+'Points/locations.pkl')
sample_locations = cPickle.load(fname)
fname.close()

fname = open('./Data/'+str(nPoints)+'Points/interpolatedData.pkl')
data_interp = cPickle.load(fname)
fname.close()

plt.ion()
plt.figure(1)
for time in range(data.shape[0]):
    plt.xlabel('Time Slice ' + str(time) + '/'+str(data.shape[0]))
    
    
    plt.subplot(3,1,1)
    plt.ylabel('Model Data')
    plt.imshow(data[time,:,:],vmin=-11.5,vmax=30.5)
    plt.axis([0, 180 ,0 ,45])
    
    plt.subplot(3,1,2)
    plt.ylabel('Interpolated Data')
    plt.imshow(data_interp[time,:,:],vmin=-11.5,vmax=30.5)
    plt.axis([0, 180 ,0 ,45])
    
    plt.subplot(3,1,3)
    plt.ylabel('Interp Data + Loc')
    plt.imshow(data_interp[time,:,:],vmin=-11.5,vmax=30.5)
    
    ind = where(sample_locations[:,0] == time)[0]
    plt.plot(sample_locations[ind,2],sample_locations[ind,1],'k.')
    plt.axis([0, 180 ,0 ,45])
    
    
    plt.draw()

    
