#variance_seasonal_cycle.py

import numpy as np
from numpy import ma
from scipy import stats
from scipy import ndimage
from scipy.io.netcdf import netcdf_file
import itertools

def example():
    timeseries = np.sin(np.arange(0, 20*np.pi, 0.01))
    timeseries = timeseries+1

    noise = np.random.uniform(0, 1, timeseries.shape)

    noisy_timeseries = timeseries+noise

    print stats.linregress(np.log(noisy_timeseries), timeseries)

nc = netcdf_file('/home/nicholas/data/netcdf_files/CFLX_2000_2009.nc', 'r')
cflux = nc.variables['Cflx'][:, :45, :180]
nc.close()

cflux = ma.masked_values(cflux, 1e20)

cflux_m = cflux
# Have not removed (masked) values above 99th percentile
# Try 
def mask_percentile(cflux_m=cflux_m):
    iterator = itertools.product(np.arange(730), np.arange(45), np.arange(180))
    for item in iterator:
        print item
        if cflux_m[item] > stats.scoreatpercentile(cflux_m[:, item[1], item[2]], 99):
            cflux_m[item] = ma.masked
    return cflux_m
               
mean_seasonal_cycle = ma.mean(np.split(cflux_m, 10, axis=0), axis=0)
mean_seasonal_cycle_repeat = np.tile(mean_seasonal_cycle, (10, 1, 1))
mscr_filtered = ndimage.filters.gaussian_filter1d(mean_seasonal_cycle_repeat, 1, axis=0)


