# -*- coding: utf-8 -*-
import cPickle

from numpy import *
from random import *
from scipy import sparse,linalg
from scipy.sparse import linalg # changed by NP
import os, sys


fname = open('data_co2.pkl', 'r')
data = cPickle.load(fname)
#t1 = data[0,:,:]


from time import *

valid_locations = []
tot_valid = 0
for t in range(data.shape[0]):
    for lon in range(data.shape[1]):
        for lat in range(data.shape[2]):
            if abs(data[t,lon,lat])>1e-10:
                tot_valid +=1
                
locations = zeros((tot_valid,3))
values = zeros((tot_valid))           



counter = 0
for t in range(data.shape[0]):
    for lon in range(data.shape[1]):
        for lat in range(data.shape[2]):
            if abs(data[t,lon,lat])>1e-10:
                values[counter] = data[t,lon,lat]
                locations[counter,0] = t
                locations[counter,1] = lon
                locations[counter,2] = lat
                valid_locations.append(counter)
                counter +=1
                




temp_valid_locations = valid_locations[:]
errors = []

all_pos = []
n_sample_points = 4
sample_locations = zeros((n_sample_points,3))
sample_values = zeros((n_sample_points))
for i in range(n_sample_points):
    ind =  int(random()*len(temp_valid_locations))
    sample_locations[i,:] = locations[temp_valid_locations[ind],:]
    sample_values[i] = values[temp_valid_locations[ind]]
    temp_valid_locations.pop(ind)
    all_pos.append(ind)
    


from FModule import *
for i in range(7400):
    print "=============================="
    print "iteration ",i

    r = 20.0
    
    r_vector,mat = fortranmodule.rbf_function_matrix(r,sample_values,sample_locations,locations)

    r1 = sparse.linalg.spsolve(mat.copy(),r_vector.copy())

    
    s_out = fortranmodule.rbf_evaluate(r,sample_values,sample_locations,locations,r1)[:,0]
  
    
    

    errors = abs(s_out - values)
    max_error = max(errors)
    tot_error = sum((s_out - values)**2)
    pos = where(errors==max_error)[0]




    n_sample_points += 1
    sample_locations_temp = sample_locations.copy()
    sample_values_temp = sample_values.copy()

    sample_locations = zeros((n_sample_points,3))
    sample_values = zeros((n_sample_points))

    sample_values[0:n_sample_points-1] = sample_values_temp
    sample_locations[0:n_sample_points-1,:] = sample_locations_temp

    sample_values[n_sample_points-1] = values[pos]
    sample_locations[n_sample_points-1,:] = locations[pos,:]


    execfile('output.py')
    
 
