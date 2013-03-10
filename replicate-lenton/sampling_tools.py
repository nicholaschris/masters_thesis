#!usr/bin/python
#sampling_tools.py
'''
NOTES:
Copy and paste from /home/nicholas/projescts/skeleton/replicate_lenton/sampling_tools.py
'''

import numpy as np
from numpy import ma 
from scipy.io.netcdf import netcdf_file 
import itertools
import timeit
from matplotlib import pyplot as plt
import os
import cPickle

### ==============================================================================================================
'''
Use this if using seperately from the lenton class.
'''
def get_an_array(nc_var):
    # input directory 
    indir = os.getcwd() + '/../data/netcdf_files/'
    ### get the flux data
    print "Fetching data from netcdf file..."
    #nc = netcdf_file(indir + 'CFLX_2000_2009.nc','r')
    nc = netcdf_file(indir + 'cflx_dpco2_1998_2007.nc','r')
    data  = nc.variables['Cflx'][:,0:45,0:180]
    #~ cflux  = ma.masked_values(cflux,1e+20)
    nc.close()
    print "Fetching data from netcdf file... Done"
    print "Fetching mask from netcdf file..."
    nc_mask=netcdf_file(indir + 'ORCA2_landmask.nc','r')
    mask = ma.masked_values(nc_mask.variables['MASK'][:, 0, :45, :180], -1e34)
    nc_mask.close()
    print "Fetching mask from netcdf file... Done"
    print "Fetching grid from netcdf file..."
    nc_grid=netcdf_file(indir + 'ORCA2.0_grid.nc','r')
    grid_area = ma.masked_values(nc_grid.variables['area'][:45, :180], -1e34)
    grid_area = grid_area*mask
    grid_lat = ma.masked_values(nc_grid.variables['lat'][:45, :180], -1e34)
    grid_lon = ma.masked_values(nc_grid.variables['lon'][:45, :180], -1e34)
    nc_mask.close()
    print "Fetching grid from netcdf file... Done"
    data = data*mask
    '''
    #data = data*1e09
    cflux = data
    cflux = cflux/(1e15) # petagrams36
    cflux = cflux*(31536000) # seconds in year
    cflux = cflux*5.9076159e+13 # total area of Southern Ocean
    cflux = cflux*12 # mol
    data = cflux
    '''
    one_year_array = data[:360, :, :]
    time_steps, lat_steps, lon_steps = np.shape(one_year_array)
    year_stack = np.split(data, 10, axis=0)
    year_stack = ma.array(year_stack)
    #year_stack = ma.masked_where(year_stack > 1e+10, year_stack)
    annual_cycle = ma.mean(year_stack, axis=0)[:360, :, :]
    return annual_cycle, data, year_stack, grid_area, grid_lat, grid_lon

annual_cycle, cflux, year_stack, grid_area, grid_lat, grid_lon = get_an_array('Cflx')

###==============================================================================

area = np.sum(grid_area)
print "Area of Southern Ocean is ", area
unit_conversion = 1000 * 86400 # 12 / 1e15 * area * 31536000  # mol * area * seconds / petagrams

def sample_all_realizations(year_array, time_frequency, lon_frequency):
	'''
	Given a sampling frequency, sample all the grid realizations at that frequency.
	'''
	print "Sampling all the realisations..."
	time_steps, lat_steps, lon_steps = np.shape(year_array)
	lat_indices = np.arange(0, lat_steps, 1)
	list_of_means = []
	list_of_stdevs = []
	time_indices = np.arange(0, time_steps, time_frequency)
	for time_count in np.arange(time_frequency):
		lon_indices = np.arange(0, lon_steps, lon_frequency)
		for lon_count in np.arange(lon_frequency):
			grid = ma.take(ma.take(ma.take(year_array, \
						(lon_indices), axis=2), \
						(lat_indices), axis=1), \
						(time_indices),axis=0)
			list_of_means = np.append(list_of_means, np.mean(grid)*unit_conversion)
			lon_indices = lon_indices - 1
		time_indices = time_indices - 1
	return ma.mean(list_of_means), ma.std(list_of_means)

def sampler_fast(array):
	'''
	Sample the grid possibilities but skip freq which results in same number of sample tracks.
	'''
	print "Forcing all the strategies..."
	time_steps, lat_steps, lon_steps = np.shape(array)
	Ax = [time_steps/i for i in range(1, time_steps+1, 1) if time_steps%i == 0]
	print Ax
	Cx = [lon_steps/i for i in range(1, lon_steps+1, 1) if lon_steps%i == 0]
	print Cx
	M = np.zeros((np.size(Ax), np.size(Cx)))
	a_ind = 0
	for a in Ax:
		c_ind = 0
		for c in Cx:
			M[a_ind, c_ind] = (sample_all_realizations(array, int(a), int(c))[1])*(2/np.sqrt(10))
			c_ind += 1
		a_ind += 1
	return M, Cx, Ax
	
def sampler_slow(array):
	'''
	Sample all the grid possibilities.
	'''
	return M	

if __name__ == '__main__':
	#~ sample = 'sample_all_the_years_at_90_15'
	#~ sample_all_the_years_at_90_15 = [sample_all_realizations(year_stack[x, :360, :, :], 90, 15)[1] for x in range(10)]
	sam = np.round(sampler_fast(annual_cycle), 2)
	sample_errors = np.round(sam[0], 2)
	Cx = sam[1]*2
	Ax = sam[2]
	#~ np.savetxt(os.getcwd() + '/output/sample_errors-'+sample+'.csv', sample_errors, delimiter=',')
	f = open(os.getcwd() + '/output/sample_errors.pkl', 'wb')
	cPickle.dump(sample_errors, f)
	f.close()
	plt.close('all')
	plt.ylabel('$Sampling$ $frequency$ $(days)$')
	plt.xlabel('$Sampling$ $frequency$ $(^{\circ}$ $lat)$')
	plt.xticks( np.arange(np.size(Cx)), Cx, rotation=45)
	plt.yticks( np.arange(np.size(Ax)), Ax)
	plt.pcolormesh(sample_errors, vmax=0.5); plt.axis('tight'); plt.colorbar()
	name='cflx'
	plt.savefig(os.getcwd() + '/output/new_brute_pcolormesh_'+name+'.png')
	plt.close('all')
	plt.ylabel('$Sampling$ $frequency$ $(days)$')
	plt.xlabel('$Sampling$ $frequency$ $(^{\circ}$ $lat)$')
	plt.xticks( np.arange(np.size(Cx)), Cx, rotation=45)
	plt.yticks( np.arange(np.size(Ax)), Ax)
	plt.contourf(sample_errors, np.arange(0, 0.52, 0.02)); plt.axis('tight'); plt.colorbar()
	name='cflx'
	plt.savefig(os.getcwd() + '/output/new_brute_countourf_'+name+'.png')
	




