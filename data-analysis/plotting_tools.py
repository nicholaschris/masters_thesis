#!usr/bin/python

#plotting_tools.py

'''
The plotting tools for the replicate_lenton application.
'''

import os
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import ma
from scipy.io.netcdf import netcdf_file
import pylab

from mpl_toolkits.basemap import Basemap
from matplotlib import cm

import scipy as sp
import datetime

import inspect, os
import sys

pathname = os.path.dirname(sys.argv[0])  

#~ plt.xticks(horizontalalignment='left')
#~ plt.yticks(verticalalignment='bottom')


mpl.rcParams.update({'font.size': 8})
#~ mpl.rcParams['xtick.direction'] = 'out'
#~ mpl.rcParams['ytick.direction'] = 'out'

FIGDIR = '~/home/Desktop/' #/home/nicholas/masters/figures/new_lenton_plots/'
NC_PATH = '/home/nicholas/data/netcdf_files/'
#~ NC_FILE =  'CFLX_2000_2009.nc'
NC_FILE = 'cflx_dpco2_1998_2007.nc'

unit_changer = -(365*24*60*60) * 12 / 1e12
longname = 'CFLX_2000_2009' #'cflx_dpco2_1998_2007'

NC_VAR_1 = 'Cflx'
longname_1 = 'DICflux'
NC_VAR_2= 'Dpco2'
longname_2 = 'dpco2'

mask_value = 1e20

### Get netcdf variable data
nc = netcdf_file(NC_PATH + NC_FILE,'r')
data_cflux = ma.masked_values(nc.variables[NC_VAR_1][:, :50, :], mask_value) # Can be till 40 or 50
nc.close()

nc_grid = netcdf_file('/home/nicholas/thesis/data/netcdf_files/ORCA2.0_grid.nc', 'r')
area = nc_grid.variables['area'][:50, :]
nc_grid.close()

time_end = np.shape(data_cflux)[0]

data_cflux_new = np.zeros(np.shape(data_cflux))

for time in np.arange(time_end):
	data_cflux_new[time, :, :] = data_cflux[time, :, :] * area
	
data_cflux_new = data_cflux_new * unit_changer



### Get netcdf variable data
nc = netcdf_file(NC_PATH + NC_FILE,'r')
data_dpco2 = ma.masked_values(nc.variables[NC_VAR_2][:, :50, :], mask_value)
nc.close()


### get landmask
nc = netcdf_file('/home/nicholas/data/netcdf_files/ORCA2_landmask.nc','r')
mask = ma.masked_values(nc.variables['MASK'][:, :, :50, :], -9.99999979e+33)
nc.close()

#~ vmin=-75
#~ vmax=75

#~ vmin_dpco2  = -75 #np.min(data)
#~ vmax_dpco2 = 75 #np.max(data)

#~ ### Mask the dat arrays
#~ data_cflux = data_cflux*mask[0, :, :, :]
data_dpco2 = data_dpco2*mask[0, :, :, :]

#~ if longname_1 == 'DICflux':
	#~ data_cflux = data_cflux*1e09 # convert to nanomol
	#~ vmin_cflux  = -100 #np.min(data)
	#~ vmax_cflux = 100 #np.max(data)
	
### Split cflux data into years	
year_stack_cflux = np.split(data_cflux, 10, axis=0)
year_stack_cflux = ma.array(year_stack_cflux)
#~ print np.shape(year_stack)

#~ ### Split dpco2 data into years
#~ year_stack_dpco2 = np.split(data_dpco2, 10, axis=0)
#~ year_stack_dpco2 = ma.array(year_stack_dpco2)

### ---------------------###
### REGRID DATA ###
### ---------------------###
nc_grid = netcdf_file(NC_PATH+ 'ORCA2.0_grid.nc','r')
lon = nc_grid.variables['lon'][0:50,:]
lat = nc_grid.variables['lat'][0:50,:]
area = nc_grid.variables['area'][0:50,:]
mask = nc_grid.variables['mask'][0,0:50,:]
nc_grid.close()

lon_min = lon.copy()
i,j = np.where(lon_min >= 180.) # elements of lon_min that are over 180
lon_min[i,j] = lon_min[i,j] - 360. # takes those elements and subtracts 360 from them

### ==============================================================================================================
### get rid of the funny extra lon and do the same for the lat array ! 
iw = np.where(lon_min[0,:] >= lon_min[0][0])[0] # are the elements that are greater or equal to the first element ie. 78.000038
ie = np.where(lon_min[0,:] < lon_min[0][0])[0] # are the elements less than 78.000038

### puts the lon in order from -180 to 180 and removes the extra 80 at the end
lon = np.concatenate((np.take(lon_min,ie,axis=1),np.take(lon_min,iw,axis=1)),axis=1)[:,:-1]
lat = np.concatenate((np.take(lat,ie,axis=1),np.take(lat,iw,axis=1)),axis=1)[:,:-1]

# The data that is to be plotted needs to be regridded
### Regrid cflux data
bm_array_cflux = [ma.concatenate((ma.take(data_cflux_new[i, :, :],ie,axis=1),ma.take(data_cflux_new[i, :, :],iw,axis=1)),axis=1)[:,:-1] for i in range(time_end)]
bm_array_cflux = ma.array(bm_array_cflux)
#~ return bm_array
#~ self.regridded_array = ma.masked_values(bm_array, 1e+20)

#~ ### Regrid dpco2 data
#~ bm_array_dpco2 = [ma.concatenate((ma.take(data_dpco2[i, :, :],ie,axis=1),ma.take(data_dpco2[i, :, :],iw,axis=1)),axis=1)[:,:-1] for i in range(3650)]
#~ bm_array_dpco2 = ma.array(bm_array_dpco2)

years = np.arange(1998, 2008, 1)


year_stack_bmcflux = np.split(bm_array_cflux, 10, axis=0)
year_stack_bmcflux = ma.array(year_stack_bmcflux)

#~ year_stack_bmdpco2 = np.split(bm_array_dpco2, 10, axis=0)
#~ year_stack_bmdpco2 = ma.array(year_stack_bmdpco2)

ll_lon = np.floor(np.min(lon)) 
ll_lat = np.floor(np.min(lat)) 
ur_lon = np.ceil(np.max(lon))
ur_lat = np.ceil(np.max(lat))

lon_list = np.round(lon[0])

def plot_decadal_mean(array, lat=lat, lon=lon):
	'''
	Plotting function to plot the decadel mean of the model data.
	'''
	vmin = -1.5
	vmax = 1.5
	plt.close('all')
	
	cmap = cm.RdBu_r
	
	#~ plt.figure(figsize=(15, 10), dpi=300)
	
	lon_0 = 180.
	lat_0 = -90.
	bounding_lat = np.ceil(np.max(lat))

	curr_map = Basemap(projection='spstere',boundinglat=bounding_lat,lon_0=lon_0,resolution='h', area_thresh=50.)
	x, y = curr_map(lon, lat)
	
	decadal_mean = ma.mean(array, axis=0)
	
	im = curr_map.pcolormesh(x, y, decadal_mean, vmin=vmin, vmax=vmax, cmap=cmap); plt.axis('tight'); # plt.colorbar()

	#draw parallels and meridians.
	delat = 10.
	circles = np.arange(-90.,-20., delat)
	curr_map.drawparallels(circles, labels=[0,0,0,0], fontsize=12, linewidth=1.5, color='k')
	
	delon = 20.
	meridians = np.arange(-180.,180., delon)
	curr_map.drawmeridians(meridians, labels=[1,1,0,1], fontsize=12, linewidth=1, color='k')
	
	#draw coastlines and continents
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')

	cb = plt.colorbar(im, orientation = 'horizontal', aspect = 55, pad = 0.08, drawedges=False)
	for t in cb.ax.get_yticklabels():
		t.set_fontsize(12)
		t.set_color('k')
		
	# put file path and name on figure
	plt.savefig('/home/nicholas/thesis/data-analysis/cfluxdailydecmean.png')



#~ multiplot_average_axes()

#~ plot_2dft()

#~ plot_ft_timeseries()

#~ plot_timeseries()

#~ plot_time_lon()

#~ plot_time_lat()

#~ plot_stdev_interannual()
#~ multiplot_stdev_interannual()

#~ plot_stdev_all_data()
#~ multiplot_stdev_all_data()

#~ plot_stdev_seasonal_cycle()
#~ multiplot_stdev_seasonal_cycle()

#~ plot_stdev_non_seasonal()
#~ multiplot_stdev_non_seasonal_cycle()

#~ plot_summer_winter_means()
	
#~ plot_annual_means()

plot_decadal_mean(bm_array_cflux)
#~ multiplot_decadal_mean()


def regrid_array(data=data_cflux):
	'''
	#Could be put with plotting tools???
	# Regrid array to be used with Basemap
	# Only works if the same latitudes and longitudes are selected from netdcf file and grid
	# Uses the ORCA netcdf file
	### transform the longitude of ORCA onto something that basemap can read
	### The ORCA grid starts at 80 and goes to 440
	### What we want: starts at 80 and goes to 180 and then switches to -180 and goes to 80
	### this method 
	'''
	from Scientific.IO.NetCDF import NetCDFFile
	#nc_grid_file = choose_netcdf_file()
	#~ indir = raw_input('Where is the ORCA netcdf file located? \n')
	nc_grid = NetCDFFile(NC_PATH+ 'ORCA2.0_grid.nc','r')
	lon = nc_grid.variables['lon'][0:40,:]
	lat = nc_grid.variables['lat'][0:40,:]
	area = nc_grid.variables['area'][0:40,:]
	mask = nc_grid.variables['mask'][0,0:40,:]
	nc_grid.close()
	
	lon_min = lon.copy()
	i,j = np.where(lon_min >= 180.) # elements of lon_min that are over 180
	lon_min[i,j] = lon_min[i,j] - 360. # takes those elements and subtracts 360 from them

	### ==============================================================================================================
	### get rid of the funny extra lon and do the same for the lat array ! 
	iw = np.where(lon_min[0,:] >= lon_min[0][0])[0] # are the elements that are greater or equal to the first element ie. 78.000038
	ie = np.where(lon_min[0,:] < lon_min[0][0])[0] # are the elements less than 78.000038

	### puts the lon in order from -180 to 180 and removes the extra 80 at the end
	lon = np.concatenate((np.take(lon_min,ie,axis=1),np.take(lon_min,iw,axis=1)),axis=1)[:,:-1]
	lat = np.concatenate((np.take(lat,ie,axis=1),np.take(lat,iw,axis=1)),axis=1)[:,:-1]

	# The data that is to be plotted needs to be regridded
	bm_array = [ma.concatenate((ma.take(data[i, :, :],ie,axis=1),ma.take(data[i, :, :],iw,axis=1)),axis=1)[:,:-1] for i in range(3650)]
	bm_array = ma.array(bm_array)
	return bm_array
	#~ self.regridded_array = ma.masked_values(bm_array, 1e+20)


#~ bm_data = regrid_array()

def create_daily_calendar(origin='days since 1998-01-01 00:00:00',cal_length=3652):
	'''
	create a daily calendar given origin and number of days ... 
	'''
	import numpy as np
	from mpl_toolkits.basemap import num2date	
	time = np.arange(0,cal_length)
	origin =  origin
	fdates = num2date(time[:],origin)
	years  = np.asarray([fdate.year for fdate in fdates])[...,np.newaxis]
	months  = np.asarray([fdate.month for fdate in fdates])[...,np.newaxis]
	days  = np.asarray([fdate.day for fdate in fdates])[...,np.newaxis]
	daily_calendar = np.concatenate((years,months,days),axis=1)
	return daily_calendar



