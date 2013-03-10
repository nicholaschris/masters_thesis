#!usr/bin/python
# ~/Desktop/ga-session-2012-11-12/load_cflux_masked.py

'''This module loads the cflux netcdf file
 or any netcdf file that needs to be masked.

'''

import numpy as np
from numpy import ma
from scipy.io.netcdf import netcdf_file
import os

cflux_5day = ''
cflux_daily = ''
dpco2_daily = ''
pco2_daily = ''


cflux_5day = '/data/netcdf_files/CFLX_2000_2009.nc'
cflux_dpco2_daily = '/data/netcdf_files/cflx_dpco2_1998_2007.nc'
### default values (These can be altered in make_gene_map.py)
file_name = os.getcwd() + '/../' + cflux_5day # cflux_dpco2_daily
time_start = 0
time_end = 730 #730 for 5day #3650 for daily
lat_start = 0
lat_end = 40 # for southern Ocean
lon_start = 0
lon_end = 182
masked_value = 1e+20

unit_changer = 60*60*24 * 1000 # micromol carbon per m per day (86400000)


def load_file(file_name = file_name, time_start = time_start, 
		      time_end = time_end, lat_start = lat_start, lat_end = lat_end,
	              lon_start = lon_start, lon_end = lon_end, masked_value = masked_value):
	nc = netcdf_file(file_name, 'r')
	new_array = nc.variables['Cflx'][time_start:time_end, lat_start:lat_end, lon_start:lon_end]
	nc.close()
	new_array = ma.masked_values(new_array, masked_value)
	new_array = ma.array(new_array, dtype=np.float32)
	new_array = new_array * unit_changer
	return new_array    

if __name__ == '__main__':
	array = load_file()

