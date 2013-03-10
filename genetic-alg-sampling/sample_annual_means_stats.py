#! usr/bin/python
# ~/Desktop/ga-session-2012-11-12/sample_annual_means.py

#TO DO 
# Fix snapshot of data, why doesn't it get correct values??? FIXED


import sys
import cPickle
import numpy as np
from numpy import ma
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
import os

from scipy.io.netcdf import netcdf_file

import datetime
### Some stuff to save time and date on plot and filenames
dt = datetime.datetime.now()
date_time_str = dt.ctime()
date_time_short = date_time_str[11:]
print date_time_short

# Set parameters
time_end = 730
lat_end = 45
lon_end = 180

unit_changer = 60*60*24 * 1000 # micromol carbon per m per day (86400000)

output_dir = os.getcwd() + '/output/plots/sample_annual_means/'
coord_dir = "/output/coords/"
coords = "coords-DecadalMean-2013-01-14_18:45-49-100000-stats.pkl"
coord_filename = os.getcwd() + coord_dir + coords
#coord_filename = os.getcwd() + "/output/coords/Coordinates_Values_2012-12-17-8h56.pkl"
#coord_filename = sys.argv[1]
#~ coord_filename = raw_input("Which cPickle file do you want to open?: \n")
f = open(coord_filename, 'rb')
coords_and_values = cPickle.load(f)
f.close()

data_cflux_5day_name = "data_cflux_5day" #"data_CFLX_5day_ave"
'''
data_filename = os.getcwd() + '/../data/pkl_files/' + data_cflux_5day_name + '.pkl' 
#data_filename = sys.argv[2]
#~ data_filename = raw_input("Which cPickle file do you want to open?: \n")
#~ nc = netcdf_file()
f = open(data_filename, 'rb')
data_cflux_5day = cPickle.load(f)
f.close()

data_cflux_daily_name = "data_cflux" #"data_cflx_daily_unmasked"
data_filename = os.getcwd() + '/../data/pkl_files/' + data_cflux_daily_name + '.pkl' 
#data_filename = sys.argv[2]
#~ data_filename = raw_input("Which cPickle file do you want to open?: \n")
#~ nc = netcdf_file()
f = open(data_filename, 'rb')
data_cflux_daily = cPickle.load(f)
f.close()

data_dpco2_daily_name = "data_dpco2" #"data_dpco2_daily_unmasked"
data_filename = os.getcwd() + '/../data/pkl_files/' + data_dpco2_daily_name + '.pkl' #sys.argv[2]
#~ data_filename = raw_input("Which cPickle file do you want to open?: \n")
#~ nc = netcdf_file()
f = open(data_filename, 'rb')
data_dpco2_daily = cPickle.load(f)
f.close()
'''
masked_value = 1e20
nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/CFLX_2000_2009.nc', 'r')
data_cflux_5day = nc.variables['Cflx'][:, :lat_end, :lon_end]
nc.close()
data_cflux_5day = ma.masked_values(data_cflux_5day, masked_value)
data_cflux_5day = ma.array(data_cflux_5day, dtype=np.float32)



data_name = data_cflux_5day_name
data = data_cflux_5day*unit_changer
year_stack = np.split(data, 10)
year_stack = ma.array(year_stack)
print "Year stack has shape: ", np.shape(year_stack)

decadal_mean = ma.mean(data, 0)
dec_mean = ma.mean(decadal_mean)
dec_stdev = ma.std(decadal_mean)
dec_range = np.abs(ma.max(decadal_mean) - ma.min(decadal_mean))
samples = []
for item in coords_and_values:
    samples.append(decadal_mean[ item[1], item[2]])

samples = ma.array(samples)
samples_mean = ma.mean(samples)
samples_stdev = ma.std(samples)
samples_range = np.abs(ma.max(samples) - ma.min(samples))
original_fitness = np.abs(dec_mean - samples_mean) + np.abs(dec_stdev - samples_stdev) \
 + np.abs(dec_range - samples_range) 
   

#~ x = 0
year_sample_dict_data = {}
year_sample_list_data = {}
for x in np.arange(np.shape(year_stack)[0]):
	year_value_list = []
	year_values_dict = {}
	x_nodes = []
	y_nodes = []
	z_nodes = []
	for item in coords_and_values:
		x_nodes.append(item[0])
		y_nodes.append(item[1])
		z_nodes.append(item[2])
		year_values_dict[item] =  ma.mean(year_stack, axis=1)[x,item[1],item[2]]
		year_value_list = np.append(year_value_list, year_stack[x,item[0],item[1],item[2]])
	year_sample_list_data[x] = year_value_list
	year_sample_dict_data[x] = year_values_dict
	print "Creating dictionary"
	#print year_sample_dict_data[x]
	
### We have a values list and coordinates
### If using 5day averages on daily data, the time (x_nodes) is multiplied by 5


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
	


	
years = np.arange(10)

#number_of_gens = str(5)
filename = os.getcwd() + '/output/tex/fast/annual_means-'+coords[:-4]+'-' +str(np.round(original_fitness, 2))+ '.tex'
#~ filename = raw_input("Enter filename to write table to: ")
fobj = open(filename, 'w')

data_name = data_name #"CFLX5DAY"
prop_or_curr = "Proposed" # Change this depending on which sampling strategy is used...
text_start = " \
\\begin{center} \n \
\\begin{table}[h] \n \
\\centering \n \
\caption{Comparison of the total Simulated Uptake With the Uptake from the \n \
Sampling strategy using " +str(np.size(x_nodes))+ " locations with a \n \
fitness "+str(np.round(original_fitness, 2))+" and the Sampling Error Introduced for " + data_name + ".} \n \
\\begin{tabular}[tbp]{@{}llllll@{}} \n \
\\toprule \n \
\scriptsize{Year} & \
\scriptsize{\specialcell{Model\\Mean\\}} \n \
& \scriptsize{\specialcell{Model\\Standard\\Deviation}} \n \
& \scriptsize{\specialcell{Sample\\Mean}} \n \
& \scriptsize{\specialcell{Sample\\Standard\\Deviation}} \n \
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

