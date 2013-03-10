# process_locations.py

import os
import cPickle
import numpy as np
from scipy.interpolate import Rbf
from scipy.io import netcdf
from numpy import ma
from matplotlib import pyplot as plt

def use_netcdf_files():
	nc = netcdf.netcdf_file('/home/nicholas/data/netcdf_files/CFLX_2000_2009.nc', 'r')
	all_data = nc.variables['Cflx'][:, :45, :180]
	nc.close()
	all_data = all_data * 1000 * 24 * 60 * 60
	all_data = ma.masked_values(all_data, 1e20)

	nc = netcdf.netcdf_file('/home/nicholas/data/netcdf_files/ORCA2.0_grid.nc', 'r')
	mask = nc.variables['mask'][0, :45, :180]
	nc.close()
	mask = ma.masked_values(mask, -1e34)
	return all_data, mask

dir_name = '7000Points'
try:
    os.mkdir(os.getcwd() + '/process_loc_plots/' + dir_name)
except OSError:
    print "Dir exists"
file_name = os.getcwd() + '/Data/' + dir_name + '/locations.pkl'
with open(file_name, 'r') as f:
	locations = cPickle.load(f)
f.closed

data_file = 'all_data_cflux5.pkl'
with open(data_file, 'r') as f:
    all_data = cPickle.load(f)
f.closed
all_data = all_data * 1000 * 24 * 60 * 60 #mmol per m**2 per day

mask_file = 'landmask.pkl'
with open(mask_file, 'r') as f:
    mask = cPickle.load(f)
f.closed
mask = ma.masked_values(mask, -1e34)

all_data = all_data * mask

data_split = np.split(all_data, 10)
data_split = np.array(data_split)
data_split = ma.masked_values(data_split, 1e20)
for year in range(10):
    for t in range(73):
        data_split[year, t, :, :] = data_split[year, t, :, :] * mask

values = {}
lat_lon_loc = {}

# Locations data from Main.py does contain some "masked locations". 
# Needs to be adressed in Main.py

for t in range(73):
    lat_lon_loc[t] = []
    for item in locations:
        if ma.is_masked(data_split[0, t, item[1], item[2]]):
            print "MASKED"
        else:
            if item[0] == t:
                lat_lon_loc[t].append([item[1], item[2]])            
    lat_lon_loc[t] = np.array(lat_lon_loc[t])
    
for year in range(10):
    for t in range(73):
        values[(year, t)] = []
        for item in lat_lon_loc[t]:
            if ma.is_masked(data_split[year, t, item[0], item[1]]):
                print "MASKED"
            #~ else:
                #~ if item[0] == t:
            values[(year, t)].append(data_split[year, t, item[0], item[1]])            
        #~ for item in locations:
            #~ if ma.is_masked(data_split[year, t, item[1], item[2]]):
                #~ print "MASKED"
            #~ else:
                #~ if item[0] == t:
                    #~ values[(year, t)].append(data_split[year, item[0], item[1], item[2]])

# Need to get rid of masked values from values and lat_lon_loc! DONE

rbf_type="gaussian"
eps=40
interp_data = np.zeros((10, 73, 45, 180))
error = np.zeros((10, 73, 45, 180))
xx, yy = np.meshgrid(np.arange(45), np.arange(180))
xi = range(45)
yi = range(180)
for year in range(10):
    for t in range(73):
        print (year, t)
        #~ print lat_lon_loc[t][:, 0], lat_lon_loc[t][:, 1], values[(year, t)]
        #try:
        rbf_func  = Rbf(lat_lon_loc[t][:, 0], lat_lon_loc[t][:, 1], values[(year, t)], function=rbf_type, epsilon=eps)
        interp_data[year, t, :, :] = np.swapaxes(rbf_func(xx, yy), 1, 0)
        interp_data[year, t, :, :] = interp_data[year, t, :, :] * mask
        error[year, t, :, :] = data_split[year, t, :, :] - interp_data[year, t, :, :]
        error[year, t, :, :] = error[year, t, :, :] * mask
        #except ValueError:
            #print "Not enough points" # WHY is it getting error: No NaNs or Inf???
            # ANSWER --> There are masked locations in values
            # BECAUSE --> The locations where not chosen correctly
            #pass
        '''
        plt.close('all')
        plt.subplot(3, 1, 1); plt.pcolormesh(data_split[year, t, :, :], vmin = -15, vmax = 15); plt.colorbar()
        plt.ylabel('Model Data')
        plt.subplot(3, 1, 2); plt.pcolormesh(interp_data[year, t, :, :], vmin = -15, vmax = 15); plt.colorbar() 
        plt.ylabel('Interpolated Data')
        plt.subplot(3, 1, 3); plt.pcolormesh(error[year, t, :, :], vmin = -15, vmax = 15); plt.colorbar()
        plt.ylabel('Interp Data + Loc')
        for item in lat_lon_loc[t]:
            plt.plot(item[1], item[0], 'k.')
        
        plt.savefig(os.getcwd() + '/process_loc_plots/' + dir_name + '/' + str(year) + '-' + str(t) + '.png')
        plt.close('all')
        '''
interp_data = interp_data*mask
error = interp_data*mask
'''
with open(os.getcwd() + '/process_loc_data/' + dir_name + 'interp.pkl', 'w') as f:
    cPickle.dump(interp_data, f)
f.closed

with open(os.getcwd() + '/process_loc_data/' + dir_name + 'error.pkl', 'w') as f:
    cPickle.dump(error, f)
f.closed
'''
interp_data_tll = np.reshape(interp_data, (730, 45, 180))
interp_timeseries = ma.mean(ma.mean(interp_data_tll, 2), 1)
data_timeseries = ma.mean(ma.mean(all_data, 2), 1)
### Plotting
plt.close('all')
plt.plot(data_timeseries, label='Original Data')
plt.plot(interp_timeseries, label='Interpolated data')
plt.axis('tight')
plt.legend()
plt.title("Original Data and Interpolated data for sampling strategy using " + dir_name)
plt.ylabel("CO2 flux(mmol per m^2 per day)")
plt.savefig(os.getcwd()+"/process_loc_data/eps"+str(eps)+str(rbf_type)+"_timeseries"+dir_name+".png")

interp_mean = ma.mean(interp_data_tll, 0)
data_mean = ma.mean(all_data, 0)
### Plotting
plt.close('all')
ax1=plt.subplot(2, 1, 1); plt.pcolormesh(data_mean, vmin = -5, vmax = 15); 
cbar1=plt.colorbar(ax=ax1); cbar1.set_label("CO2 flux(mmol per m^2 per day)")
plt.title('Original Data Mean over time axis')
ax2=plt.subplot(2, 1, 2); plt.pcolormesh(interp_mean, vmin = -5, vmax = 15); 
cbar2=plt.colorbar(ax=ax2); cbar2.set_label("CO2 flux(mmol per m^2 per day)")
plt.title('Interpolated Data Mean over time axis')
plt.savefig(os.getcwd()+"/process_loc_data/eps"+str(eps)+str(rbf_type)+"_average"+dir_name+".png")

