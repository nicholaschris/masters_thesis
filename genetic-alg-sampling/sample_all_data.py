#!usr/bin/python
# sample_all_data.py

import numpy as np
import cPickle
import os
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.io.netcdf import netcdf_file
from numpy import ma

mpl.rcParams.update({'font.size': 8})
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'

nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/CFLX_2000_2009.nc', 'r')
cflx = nc.variables['Cflx'][:, :50, :]
nc.close()
cflx = ma.masked_values(cflx, 1e20)
timeseries = ma.mean(ma.mean(cflx, 2), 1)
timeseries = timeseries*-1000*24*60*60

pkl_file = 'coords-AllData-2013-01-16_14:12-999-100000-stats'
coord_file = os.getcwd() +'/output/coords/' +pkl_file+ '.pkl'
f = open(coord_file, 'r')
coords = cPickle.load(f)
f.close()

try: 
	os.mkdir(os.getcwd() + '/output/plots/all_data/' + pkl_file)
except OSError:
	print "path exists"
	
output_dir = os.getcwd() + '/output/plots/all_data/' + pkl_file
time = []
lat = []
lon = []
values = []

for item in coords:
    time.append(item[0])
    #~ print "time ", item[0]
    lat.append(item[1])
    #~ print "lat ",  item[1]
    lon.append(item[2])
    #~ print "lon ", item[2]
    values.append(coords[item])

###    
plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.hist(time, bins=73); ax1.axis('tight')
ax1.set_xlabel('Time (Years)')
ax1.set_ylabel('Frequency of samples')
ax2 = ax1.twinx()
ax2.plot(timeseries, 'o', color='#FF2400'); ax2.axis('tight')
ax2.set_ylabel('CO2 flux (mmol/m^2/day)')
plt.xticks(np.arange(0, 730, 73), np.arange(2000, 2010, 1))
plt.savefig(output_dir +'/histogram_time.png')

plt.close('all')
plt.hist(lat, bins=40); plt.axis('tight')
plt.savefig(output_dir +'/histogram_lat.png')

plt.close('all')
plt.hist(lon, bins=180); plt.axis('tight')
plt.xticks(np.arange(0, 180, 10), np.arange(80, 440, 20))
plt.savefig(output_dir +'//histogram_lon.png')
###

###
plt.close('all')
plt.scatter(lon, lat, 10, values); plt.colorbar(); plt.axis('tight')

plt.savefig(output_dir +'/scatter_latlon_values.png')

plt.close('all')
plt.scatter(lon, time,  10, values); plt.colorbar(); plt.axis('tight')
plt.yticks(np.arange(0, 730, 73), np.arange(2000, 2010, 1))
plt.savefig(output_dir +'/scatter_timelon_values.png')

plt.close('all')
plt.scatter(time, lat, 10, values); plt.colorbar(); plt.axis('tight')
plt.xticks(np.arange(0, 730, 73), np.arange(2000, 2010, 1))
plt.savefig(output_dir +'/scatter_timelat_values.png')
###

###
plt.close('all')
plt.scatter(lon, lat, 10, (np.array(time)%73)/6); plt.colorbar(); plt.axis('tight')
plt.savefig(output_dir +'/scatter_latlon_time.png')

plt.close('all')
plt.scatter(lon, time, 10, lat); plt.colorbar(); plt.axis('tight')
plt.yticks(np.arange(0, 730, 73), np.arange(2000, 2010, 1))
plt.savefig(output_dir +'/scatter_timelon_lat.png')

plt.close('all')
plt.scatter(time, lat, 10, lon); plt.colorbar(); plt.axis('tight')
plt.xticks(np.arange(0, 730, 73), np.arange(2000, 2010, 1))
plt.savefig(output_dir +'/scatter_timelat_lon.png')
###
