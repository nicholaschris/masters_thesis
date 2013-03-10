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
lon_freq = 30
time_start = 0
lat_start = 0
lon_start = 0
time_end = 73
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


reg_grid_samples = {}
coords = itertools.product(np.arange(time_start, time_end, time_freq), np.arange(0, lat_end, lat_freq), np.arange(0, lon_end, lon_freq))
for year in np.arange(10):
    coords = itertools.product(np.arange(time_start, time_end, time_freq), np.arange(0, lat_end, lat_freq), np.arange(0, lon_end, lon_freq))
    reg_grid_samples[year] = []
    for x in coords:
        #print year_stack[year, :, :, ][x]
        reg_grid_samples[year].append(year_stack[year, :, :, :][x])

num_loc =  np.size(reg_grid_samples[9])
print num_loc

data_mean=[]
sample_mean=[]
data_stdev=[]
sample_stdev=[]
data_range=[]
sample_range=[]
fitness = []

for year in np.arange(10):
	data_mean.append(ma.mean(year_stack[year, :, :, :]))
	sample_mean.append(ma.mean(reg_grid_samples[year]))
	data_stdev.append(ma.std(year_stack[year, :, :, :]))
	sample_stdev.append(ma.std(reg_grid_samples[year]))
	data_range.append(np.abs(ma.max(ma.mean(year_stack[year, :, :, :], 1)) - ma.min(ma.mean(year_stack[year, :, :, :], 1))))
	sample_range.append(np.abs(ma.max(reg_grid_samples[year]) - ma.min(reg_grid_samples[year])))
	fitness.append(np.abs(data_mean[year]-sample_mean[year]) + \
	np.abs(data_stdev[year]-sample_stdev[year]) + \
	np.abs(data_range[year]-sample_range[year]))
	
mean_data = data_mean
stdev_data = data_stdev
mean_interp_data = sample_mean
stdev_interp_data = sample_stdev


seasonal_cycle = ma.mean(np.split(data, 10, axis=0), axis=0)
dec_mean = ma.mean(seasonal_cycle)
dec_stdev = ma.std(seasonal_cycle)
dec_range = np.abs(ma.max(seasonal_cycle) - ma.min(seasonal_cycle))
samples = []
coords = itertools.product(np.arange(time_start, time_end, time_freq), np.arange(0, lat_end, lat_freq), np.arange(0, lon_end, lon_freq))
for item in coords:
    samples.append(seasonal_cycle[item])

samples = ma.array(samples)
samples_mean = ma.mean(samples)
samples_stdev = ma.std(samples)
samples_range = np.abs(ma.max(samples) - ma.min(samples))
original_fitness = np.abs(dec_mean - samples_mean) + np.abs(dec_stdev - samples_stdev) \
 + np.abs(dec_range - samples_range) 
   
'''
def test_stats(x):
	coords_vals = year_sample_dict_data[x]
		
	x_nodes = []	
	y_nodes = []
	z_nodes = []
	values_list = []
	for item in coords_vals:
		x_nodes.append(item[0])
		y_nodes.append(item[1])
		z_nodes.append(item[2])
		values_list.append(coords_vals[item])
	xs = x_nodes
	ys = y_nodes
	zs = z_nodes
	values_list = np.array(values_list)
	_time_end = 73
	all_data = year_stack[x, :_time_end, :lat_end, :lon_end]
	### New and improved and faster!!!
	annual_mean = np.mean(all_data)
	sample_mean = np.mean(values_list)
	annual_stdev = np.std(all_data)
	sample_stdev = np.std(values_list)
	annual_max = np.max(all_data)
	annual_min = np.min(all_data)
	sample_max = np.max(all_data)
	sample_min = np.min(all_data)
	annual_range = np.abs(annual_max - annual_min)
	sample_range = np.abs(sample_max - sample_min)
	fitness = np.abs(annual_mean-sample_mean) + np.abs(annual_stdev - sample_stdev) + np.abs(annual_range - sample_range)
	return fitness, annual_mean, sample_mean, annual_stdev, sample_stdev 


list_of_fitness=[]
for x in np.arange(10):
    print x
    list_of_fitness.append(test_stats(x))

fitness = []
mean_data = []
mean_interp_data = []
stdev_data = []
stdev_interp_data = []

for item in list_of_fitness:
	fitness = np.append(fitness, item[0])
	mean_data = np.append(mean_data, item[1])
	mean_interp_data = np.append(mean_interp_data, item[2])
	stdev_data = np.append(stdev_data, item[3])
	stdev_interp_data = np.append(stdev_interp_data, item[4])
	
'''

	
years = np.arange(10)

#number_of_gens = str(5)
filename = os.getcwd() + '/output/tex/reg_grid/seasonal_cycle-'+str(np.round(original_fitness, 2))+ '_'+str(time_freq)+'_'+str(lat_freq)+'_'+str(lon_freq)+'.tex'
#~ filename = raw_input("Enter filename to write table to: ")
fobj = open(filename, 'w')

data_name = "CFLX5DAY"
prop_or_curr = "Proposed" # Change this depending on which sampling strategy is used...
text_start = " \
\\begin{center} \n \
\\begin{table}[h] \n \
\\centering \n \
\caption{Comparison of the total Simulated Uptake With the Uptake from the \n \
Sampling strategy using " +str(num_loc)+ " locations with a \n \
fitness "+str(np.round(original_fitness, 2))+" and the Sampling Error Introduced for " + data_name + ".} \n \
\\begin{tabular}[tbp]{@{}llllll@{}} \n \
\\toprule \n \
\scriptsize{Year} & \
\scriptsize{\specialcell{Model\\\Mean\\}} \n \
& \scriptsize{\specialcell{Model\\\Standard\\\Deviation}} \n \
& \scriptsize{\specialcell{Sample\\\Mean}} \n \
& \scriptsize{\specialcell{Sample\\\Standard\\\Deviation}} \n \
& \scriptsize{Fitness} \\\ \n \
\hline \n \
"

fobj.write(text_start)
'''
for year in range(10):
	fobj.write("\scriptsize{" +str(years[year] +2000)+ "} & \scriptsize{"+str(round(list_of_fitness[year][1], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][2], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][3], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][4], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][0], 2))+"} \\\\")
	fobj.write("\n")
'''
for year in range(10):
	fobj.write("\scriptsize{" +str(years[year] +2000)+ "} & \scriptsize{"+str(round(mean_data[year], 2))+"} & \scriptsize{"+str(round(stdev_data[year], 2))+"} & \scriptsize{"+str(round(mean_interp_data[year], 2))+"} & \scriptsize{"+str(round(stdev_interp_data[year], 2))+"} & \scriptsize{"+str(round(fitness[year], 2))+"} \\\\")
	fobj.write("\n")

text_end ="\hline \n \
\scriptsize{" +str(years[0] + 2000)+ "-" +str(years[9] + 2000)+ "}    &   \
\scriptsize{"+str(round(np.mean(mean_data), 2))+" $\pm$ "+str(round(np.std(mean_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(stdev_data), 2))+" $\pm$ "+str(round(np.std(stdev_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(mean_interp_data), 2))+" $\pm$ "+str(round(np.std(mean_interp_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(stdev_interp_data), 2))+" $\pm$ "+str(round(np.std(stdev_interp_data), 2))+"} & \
\scriptsize{"+str(round(np.mean(fitness), 2))+" $\pm$ "+str(round(np.std(fitness), 2))+"} \\\ \n \
\\bottomrule \n \
\end{tabular} \n \
\end{table} \n \
\end{center} \n \
"

fobj.write(text_end)
fobj.close()

print "Output written to ", filename

