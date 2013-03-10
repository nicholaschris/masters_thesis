#! usr/bin/python
# ~/Desktop/ga-session-2012-11-12/sample_seasonal_cycle.py

'''
This script takes the solution from from the genetic
algorithm that samples the seasonal cycle and tests
it on the data for each year.

It returns a fitness values for each year, as well as 
the mean value  and the standard
deviation from the sample as well as the model data.

Values get put into a csv or a tex table?

'''

#! usr/bin/python
#  sample_annual_means.py

#TO DO 
# Fix snapshot of data, why doesn't it get correct values??? FIXED

import cPickle
import numpy as np
from numpy import ma
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os
from scipy.io.netcdf import netcdf_file

unit_changer = 60*60*24 * 1000 # micromol carbon per m per day (86400000)
# Set parameters
time_len = 365
lat_end = 45
lon_end = 180


import datetime
### Some stuff to save time and date on plot and filenames
dt = datetime.datetime.now()
date_time_str = dt.ctime()
date_time_short = date_time_str[11:]
print date_time_short

data_dir = os.getcwd() + '/output/coords/'
coord_file = "coords-AnnualCycle-2013-01-15_23:10-1999-10000-stats.pkl"
coord_filename = data_dir + coord_file
#coord_filename = data_dir + sys.argv[1] # "Coordinates_Values_2012-10-30-20h56.pkl"
#~ coord_filename = raw_input("Which cPickle file do you want to open?: \n")
f = open(coord_filename, 'rb')
coords_and_values = cPickle.load(f)
f.close()

'''
data_filename = os.getcwd() + "/../data/pkl_files/data_cflux_5day.pkl"
#~ data_filename = raw_input("Which cPickle file do you want to open?: \n")
#~ nc = netcdf_file()
f = open(data_filename, 'rb')
data = cPickle.load(f)
f.close()


data_cflux_daily_name = "data_cflux" #"data_cflx_daily_unmasked"
data_filename = os.getcwd() + '/../data/pkl_files/' + data_cflux_daily_name + '.pkl' 
#data_filename = sys.argv[2]
#~ data_filename = raw_input("Which cPickle file do you want to open?: \n")
#~ nc = netcdf_file()
f = open(data_filename, 'rb')
data_cflux_daily = cPickle.load(f)
f.close()


data_name = data_cflux_daily_name
data = data_cflux_daily*1e08
data_shape = np.shape(data)
year_stack = np.split(data, 10)
year_stack = ma.array(year_stack)
print "Year stack has shape: ", np.shape(year_stack)
'''
data_cflux_5day_name = "cflux_5day"
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

'''
### get landmask
nc = netcdf_file(os.getcwd() + '/../data/netcdf_files/ORCA2_landmask.nc','r')
mask = ma.masked_values(nc.variables['MASK'][:, :time_len, :lat_end, :lon_end], -9.99999979e+33)
nc.close()
'''

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
		#print item #DEBUG
		x_nodes.append(item[0])
		y_nodes.append(item[1])
		z_nodes.append(item[2])
		year_values_dict[item] =  year_stack[x, item[0], item[1],item[2]]
		year_value_list = np.append(year_value_list, year_stack[x,item[0],item[1],item[2]])
	year_sample_list_data[x] = year_value_list
	year_sample_dict_data[x] = year_values_dict
	
	
	
### We have a values list and coordinates
def test_rbf(x):
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
	all_data = year_stack[x, :time_len, :lat_end, :lon_end]
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
	###lat = lat_end
	###lon = lon_end

	###ti_lat = np.arange(0, lat)
	###ti_lon = np.arange(0, lon)

	###lati, loni = np.meshgrid(ti_lat, ti_lon)

	###s = len(xs)




	#### Trying out 3d RBF (insert into ga_class)
	#### meshgrid 3d
	###xxx, yyy, zzz = np.lib.index_tricks.mgrid[0:time_len, 0:lat_end, 0:lon_end]


	# use RBF
	###rbf = Rbf(xs, ys, zs, values_list, function='gaussian', epsilon=4)
	#~ rbf = Rbf(zs, ys, values_list, function='gaussian', epsilon=4)
	###ZI = rbf(xxx, yyy, zzz)
	#~ ZI = rbf(loni, lati)
	#~ ZI = np.swapaxes(ZI, 1, 0)
	#~ ZI = np.reshape(ZI, (1, 40, 180))
	#~ ZI = ZI*np.tile(mask[0, :, :, :], (0, 2, 1, 1))
	###ZI = ZI*mask[0, :, :, :]
	#~ plt.close('all')

	#~ plt.subplot(1, 1, 1)
	#~ plt.pcolormesh(ma.mean(ZI, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
	#~ plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
	#~ plt.title('RBF interpolation - multiquadrics')

	#~ plt.axis('tight')
	#~ plt.colorbar()
	#~ plt.savefig('year'+str(x)+'rbf2d' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')

	#~ plt.close('all')

	###GI = year_stack[x, :time_len, :lat_end, :lon_end]
	###print np.shape(GI) # DEBUG

	###diff_array = np.sqrt((GI - ZI)**2)

	#~ plt.pcolormesh(ma.mean(diff_array, axis=0), cmap=cm.jet, vmin=0, vmax=10); 
	#~ plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
	#~ plt.axis('tight')
	#~ plt.colorbar()
	#~ plt.savefig('year'+str(x)+'diff_array' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
	#~ plt.close('all')
	###fitness = np.sqrt(np.mean((GI - ZI)**2))
	###print np.mean(GI), np.mean(ZI), np.std(GI), np.std(ZI)
	###print "The fitness of the solution is: %g" % fitness # This answer should be the same as the final fittest solution???

	#~ plt.pcolormesh(ma.mean(GI, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
	#~ plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
	#~ plt.axis('tight')
	#~ plt.colorbar()
	#~ plt.savefig('year'+str(x)+'snap_shot' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
	#~ plt.close('all')
	###return fitness, np.mean(GI), np.mean(ZI), np.std(GI), np.std(ZI)

list_of_fitness=[]
for x in np.arange(10):
	list_of_fitness.append(test_rbf(x))

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
# number_of_gens = str(5)
filename = os.getcwd() + '/output/tex/' + coord_file[:-4] + \
data_cflux_5day_name +'.tex' # + number_of_gens + 'gens.tex'
#~ filename = raw_input("Enter filename to write table to: ")
fobj = open(filename, 'w')

prop_or_curr = "Proposed" # Change this depending on which sampling strategy is used...
text_start = " \
\\begin{center} \n \
\\begin{table}[h] \n \
\\centering \n \
\caption{Comparison of the total Simulated Uptake With the Uptake from Our \n \
"+ prop_or_curr + " Sampling and the Sampling Error Introduced.} \n \
\\begin{tabular}[tbp]{@{}lp{2.5cm}p{2.5cm}p{2.5cm}p{2.5cm}p{2.5cm}@{}} \n \
\\toprule \n \
\scriptsize{Year} & \scriptsize{Model Mean} \n \
& \scriptsize{Sample Mean} \n \
& \scriptsize{Model Standard Deviation} \n \
& \scriptsize{Sample Standard Deviation} \n \
& \scriptsize{Fitness} \\\ \n \
\hline \n \
"

fobj.write(text_start)

for year in range(10):
	fobj.write("\scriptsize{" +str(years[year] +2000)+ "} & \scriptsize{"+str(round(list_of_fitness[year][1], 2))+"} & \
	\scriptsize{"+str(round(list_of_fitness[year][2], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][3], 2))+"} & \
	\scriptsize{"+str(round(list_of_fitness[year][4], 2))+"} & \scriptsize{"+str(round(list_of_fitness[year][0], 2))+"} \\\\")
	fobj.write("\n")

text_end ="\hline \n \
\scriptsize{" +str(years[0] + 2000)+ "-" +str(years[9] + 2000)+ "}    &   \
\scriptsize{"+str(round(np.mean(mean_data), 2))+" $\pm$ "+str(round(np.std(mean_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(mean_interp_data), 2))+" $\pm$ "+str(round(np.std(mean_interp_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(stdev_data), 2))+" $\pm$ "+str(round(np.std(stdev_data), 2))+"} &  \
\scriptsize{"+str(round(np.mean(stdev_interp_data), 2))+" $\pm$ "+str(round(np.std(stdev_interp_data), 2))+"} & \
\scriptsize{"+str(round(np.mean(fitness), 2))+" $\pm$ "+str(round(np.std(fitness), 2))+"} \\\ \n \
\\bottomrule \n \
\end{tabular} \n \
\end{table} \n \
\end{center} \n \
"
'''
text_end ="\hline \n \
\scriptsize{" +str(years[0] + 2000)+ "-" +str(years[9] + 2000)+ "}    \
&   \scriptsize{"+str(np.mean(mean_data ))+"} &  \
\scriptsize{"+str(np.mean(mean_interp_data))+"} & \
\scriptsize{"+str(np.mean(stdev_data))+"} &  \
\scriptsize{"+str(np.mean(stdev_interp_data))+"} \\\ \n \
\\bottomrule \n \
\end{tabular} \n \
\end{table} \n \
\end{center} \n \
"
'''

fobj.write(text_end)
fobj.close()

print "Table printed to: ", filename
### OLD STUFF???


'''
\scriptsize{" +str(years[0])+ "}    &   \scriptsize{"+str(tot_sam_upt[0])"} &  \scriptsize{"str(sam_est_upt[0])"} &  \scriptsize{"str(sam_unc[0]"} \\ \n \"
\scriptsize{" +str(years[1])+ "}    &   \scriptsize{"+str(tot_sam_upt[1])"} &  \scriptsize{"str(sam_est_upt[1])"} &  \scriptsize{"str(sam_unc[1]"} \\ \n \
\scriptsize{" +str(years[2])+ "}    &   \scriptsize{"+str(tot_sam_upt[2])"} &  \scriptsize{"str(sam_est_upt[2])"} &  \scriptsize{"str(sam_unc[2]"} \\ \n \
\scriptsize{" +str(years[3])+ "}    &   \scriptsize{"+str(tot_sam_upt[3])"} &  \scriptsize{"str(sam_est_upt[3])"} &  \scriptsize{"str(sam_unc[3]"} \\ \n \
\scriptsize{" +str(years[4])+ "}    &   \scriptsize{"+str(tot_sam_upt[4])"} &  \scriptsize{"str(sam_est_upt[4])"} &  \scriptsize{"str(sam_unc[4]"} \\ \n \
\scriptsize{" +str(years[5])+ "}    &   \scriptsize{"+str(tot_sam_upt[5])"} &  \scriptsize{"str(sam_est_upt[5])"} &  \scriptsize{"str(sam_unc[5]"} \\ \n \
\scriptsize{" +str(years[6])+ "}    &   \scriptsize{"+str(tot_sam_upt[6])"} &  \scriptsize{"str(sam_est_upt[6])"} &  \scriptsize{"str(sam_unc[6]"} \\ \n \
\scriptsize{" +str(years[7])+ "}    &   \scriptsize{"+str(tot_sam_upt[7])"} &  \scriptsize{"str(sam_est_upt[7])"} &  \scriptsize{"str(sam_unc[7]"} \\ \n \
\scriptsize{" +str(years[8])+ "}    &   \scriptsize{"+str(tot_sam_upt[8])"} &  \scriptsize{"str(sam_est_upt[8])"} &  \scriptsize{"str(sam_unc[8]"} \\ \n \
\scriptsize{" +str(years[9])+ "}    &   \scriptsize{"+str(tot_sam_upt[9])"} &  \scriptsize{"str(sam_est_upt[9])"} &  \scriptsize{"str(sam_unc[9]"} \\ \n \
'''
'''
This script takes the solution from from the genetic
algorithm that samples the decadal mean and tests
it on the annual means.

It returns a fitness values for each year, as well as 
the mean value  and the standard
deviation from the sample as well as the model data.

Values get put into a csv or a tex table?


solution

coordinates

Need to get a list of values for those coordinates
for each year.

Have a list of coordinates
get a list of values

Do This:

for item in coords_vals:
	x_nodes.append(item[0])
	y_nodes.append(item[1])
	z_nodes.append(item[2])
	values_list.append(coords_vals[item])
	year_data = np.append(year_data, year(coords_vals[item]))
	
for item in coords_vals:
	values_list.append(coords_vals[item])
	year_data = np.append(year_data, year(coords_vals[item]))
	
Then have these 4 lists

x_nodes # time
y_nodes # lat
z_nodes #lon
values_list

Need meshgrid
lati, loni = np.meshgrid(ti_lat, ti_lon) #2D
OR
xxx, yyy, zzz = np.lib.index_tricks.mgrid[0:10, 0:40, 0:180] #3D

Use Rbf

rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)

interpolated_data = rbf(xxx, yyy, zzz)

fitness = np.sqrt(np.sum((GI - ZI)**2)) where GI is original data for the year and ZI is interpolated data

So for each year:


year  	| model mean | model stdev | sample mean | sample stdev | fitness |
1998
1999
2000
    |
2007

'''
