#! usr/bin/python
# ~/Desktop/ga-session-2012-11-12/sample_regular_grid_seasonal_cycle.py

import sys
import cPickle
import numpy as np
from numpy import ma
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import itertools

from scipy.io.netcdf import netcdf_file

import datetime
### Some stuff to save time and date on plot and filenames
dt = datetime.datetime.now()
date_time_str = dt.ctime()
date_time_short = date_time_str[11:]
print date_time_short

unit_changer = 60*60*24 * 1000 # micromol carbon per m per day (86400000)

output_dir = os.getcwd() + '/output/plots/sample_annual_means/'

#SET PARAMETERS
time_freq = 72/4
lat_freq = 1
lon_freq = 15
time_start = 0
lat_start = 0
lon_start = 0
time_end = 730
lat_end = 45
lon_end = 180

masked_value = 1e20
nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/CFLX_2000_2009.nc', 'r')
data_cflux_5day = nc.variables['Cflx'][:, :lat_end, :lon_end]
nc.close()
data_cflux_5day = ma.masked_values(data_cflux_5day, masked_value)
data_cflux_5day = ma.array(data_cflux_5day, dtype=np.float32)

#data_name = data_cflux_5day_name
data = data_cflux_5day*unit_changer
year_stack = np.split(data, 10)
year_stack = ma.array(year_stack)
print "Year stack has shape: ", np.shape(year_stack)

sample = []
coords = itertools.product(np.arange(time_start, time_end, time_freq), np.arange(0, lat_end, lat_freq), np.arange(0, lon_end, lon_freq))

for x in coords:
    #print year_stack[year, :, :, ][x]
    sample.append(data[x])

sample_size = np.size(sample)

data_mean=ma.mean(data)
sample_mean=ma.mean(sample)
data_stdev=ma.std(data)
sample_stdev=ma.std(sample)
data_range=(np.abs(ma.max(data)) - ma.min(ma.mean(data)))
sample_range=(np.abs(ma.max(sample) - ma.min(sample)))
fitness=(np.abs(data_mean-sample_mean) + \
np.abs(data_stdev-sample_stdev) + \
np.abs(data_range-sample_range))
	

f = open(os.getcwd() + '/output/tex/reg_grid/sample_all_'+str(sample_size)+'.tex', 'w')

text = 'Using a regular grid with a sample size of '+str(sample_size)+', a time frequency of '+str(time_freq)+', a longitudinal frequency of '+str(2*lon_freq)+', and sampling at the latitudinal resolution, the fitness value is '+str(fitness)+'.'

f.write(text)

f.close()



