#!usr/bin/python
# -*- coding: utf-8 -*-
#data-analysis.py
"""
Created on Thu Oct 11 20:33:09 2012
Modified on Tue Dec 04 17:59 2012
@author: nicholas
NOTES:
Copy and Paste from /home/nicholas/masters_python_script.py
"""

import inspect
import os
import sys
import numpy as np
from numpy import ma
from scipy.io.netcdf import netcdf_file
from scipy import signal
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from matplotlib import pyplot as plt
import matplotlib as mpl

plt.xticks(horizontalalignment='left')
plt.yticks(verticalalignment='bottom')
mpl.rcParams.update({'font.size': 8})
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
pathname = os.path.dirname(sys.argv[0]) 
#FIGDIR = '/home/nicholas/masters/figures/newplots/'
FIGDIR = os.getcwd() + '/output/figures/'

unit_changer_1 = 24*60*60 * -1000 #mmolC/m^2/day
unit_changer_2 = -(365*24*60*60) * 12 / 1e12 #teragramC/year

def get_data():
    nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/CFLX_2000_2009.nc', 'r')
    #~ nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/cflx_dpco2_1998_2007.nc', 'r')
    cflx = nc.variables['Cflx'][:, :50, :]
    time_end = np.shape(cflx)[0]
    print time_end
    #~ dpco2 = nc.variables['Dpco2'][:, 0:50, :]
    nc.close()
    nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/ORCA2_landmask.nc','r')
    mask = ma.masked_values(nc.variables['MASK'][:, 0, :50, :], -1e34)
    nc.close()
    cflx = cflx*mask
    #~ dpco2 = dpco2*mask
    nc = netcdf_file('/home/nicholas/thesis/data/netcdf_files/ORCA2.0_grid.nc', 'r')
    area = nc.variables['area'][0:50,:]
    nc.close()
    print np.shape(cflx)
    cflx_new = np.zeros(np.shape(cflx))
    print np.shape(cflx_new)
    for time in np.arange(time_end):
        cflx_new[time, : ,:]=cflx[time, :, :]*area
        
    cflx_new = cflx_new * unit_changer_2 
    print np.shape(cflx_new)
    # Atmospheric CO2 txt file
    #~ year, atm_co2 = [], []
    #~ f = open('/home/nicholas/data/atm_co2.csv')
    #~ for i in range(10):
        #~ s = f.readline()
        #~ year.append(s[0:4])
        #~ atm_co2.append(s[6:-1])
    # get pco2 data
    #~ pco2_air = np.zeros(np.shape(dpco2))	
    #~ ys_pco2_air = np.ma.array(np.split(pco2_air, 10, axis=0))	
    #~ for x in np.arange(10):
        #~ ys_pco2_air[x, :, :, :] = atm_co2[x]
    #~ # The pCO2 data is the atmospheric CO2 minus Delta pCO2    
    #~ pco2 = np.reshape(ys_pco2_air, np.shape(dpco2)) - dpco2
    return cflx_new #, dpco2, pco2

cflx_new = get_data()

def regrid(data_array):
    nc_grid = netcdf_file('/home/nicholas/thesis/data/netcdf_files/ORCA2.0_grid.nc','r')
    lon = nc_grid.variables['lon'][0:50,:]
    lat = nc_grid.variables['lat'][0:50,:]
    area = nc_grid.variables['area'][0:50,:]
    mask = nc_grid.variables['mask'][0,0:50,:]
    nc_grid.close()
    lon_min = lon.copy()
    i,j = np.where(lon_min >= 180.)
    lon_min[i,j] = lon_min[i,j] - 360. 
    iw = np.where(lon_min[0,:] >= lon_min[0][0])[0]
    ie = np.where(lon_min[0,:] < lon_min[0][0])[0]
    lon = np.concatenate((np.take(lon_min,ie,axis=1)[:,:-1],np.take(lon_min,iw,axis=1)),axis=1)[:,:-1]
    lat = np.concatenate((np.take(lat,ie,axis=1)[:,:-1],np.take(lat,iw,axis=1)),axis=1)[:,:-1]
    lon_list = np.round(lon[0])
    lat_list = np.round(lat[:, 0])
    timesteps = np.shape(data_array)[0]
    bm_array = [ma.concatenate((ma.take(data_array[i, :, :],ie,axis=1),ma.take(data_array[i, :, :],iw,axis=1)),axis=1)[:,:-1] for i in range(timesteps)]
    bm_array = ma.array(bm_array)
    return bm_array, lon, lat
    
cflx_new, lon, lat = regrid(cflx_new)
#~ time_end = np.shape(cflx_new)[0]
#~ pco2 = regrid(pco2)

#Basemap stuff
ll_lon = np.floor(np.min(lon)) 
ll_lat = np.floor(np.min(lat)) 
ur_lon = np.ceil(np.max(lon))
ur_lat = np.ceil(np.max(lat))
cmap = cm.jet
curr_map = Basemap(projection='cyl', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='h', area_thresh=100.)
x, y = curr_map(lon, lat)

vmin=-1.5 #250
vmax=1.5 #450

# make latitude and longitude lists
_e2w = np.arange(-180, 181, 1)
e2w= np.append(['$'+str(abs(i))+'^{\circ}E$' for i in _e2w if i<=0],['$'+str(abs(i))+'^{\circ}W$' for i in _e2w if i>0] )
_s2n = np.arange(-90, 91, 1)
s2n = np.append(['$'+str(abs(i))+'^{\circ}S$' for i in _s2n if i<=0],['$'+str(abs(i))+'^{\circ}N$' for i in _s2n if i>0] )

lon_list = np.round(lon[0])
lat_list = np.round(lat[:, 0])

#plt.yticks(np.arange(-80, -20, 10), lat_list.take(np.arange(0, 50, 10)))
#plt.xticks( np.arange(-180, 180, 60), lon_list.take(np.arange(0, 180, 30)))

def plot_tll(data, name):
    #plot decadel mean with lat/lon aves cflx dpco2
    #~ plt.figure(figsize=(10, 10), dpi=50)
    plt.subplot(221)
    im = curr_map.pcolormesh(x, y, ma.mean(data, 0), vmin=vmin, vmax=vmax, cmap=cmap); plt.axis('tight');
#    #draw parallels and meridians.
#    delat = 10.
#    circles = np.arange(-90.,-20., delat)
#    curr_map.drawparallels(circles, labels=[1,0,0,0], fontsize=12, linewidth=1.5, color='k')
#    delon = 60.
#    meridians = np.arange(-180.,180., delon)
#    curr_map.drawmeridians(meridians, labels=[1,0,0,1], fontsize=12, linewidth=1, color='k')    
    curr_map.drawcoastlines()
    curr_map.fillcontinents(color='grey',lake_color='aqua')
    #~ plt.xlabel('$Longitude$')
    #~ plt.ylabel('$Latitude$')
    plt.xticks( np.arange(-160, 180, 60), e2w[20:-1:60])
    #~ plt.xticks( np.arange(-160, 180, 60), lon_list.take(np.arange(9, 180, 30)))
    plt.yticks([-78+16, -78+16+11, -78+16+11+9, -78+16+11+9+7], s2n[20:70:10])
    #~ plt.yticks([-78+16, -78+16+11, -78+16+11+9, -78+16+11+9+7], [lat_list[16], lat_list[16+11], lat_list[16+11+9], lat_list[16+11+9+7]])
    plt.subplot(223)
    plt.pcolormesh(ma.mean(data, 1), vmin=vmin, vmax=vmax, cmap=cmap); plt.axis('tight')
    #~ plt.xlabel('$Longitude$')
    #~ plt.ylabel('$Time$')
    plt.xticks( np.arange(9, 180, 30), e2w[20:-1:60])
    plt.yticks(np.arange(0, 3650, 365), np.arange(1998, 2008), rotation=45)
    plt.subplot(222)
    plt.pcolormesh(np.swapaxes(ma.mean(data, 2), 1, 0), vmin=vmin, vmax=vmax, cmap=cmap); plt.axis('tight')
    #~ plt.xlabel('$Time$')
    #~ plt.ylabel('$Latitude$')
    plt.xticks(np.arange(0, 3650, 365), np.arange(1998, 2008), rotation=-45)
    plt.yticks([16, 16+11, 16+11+9, 16+11+9+7], s2n[20:70:10])
    plt.subplot((224), visible=False)
    plt.colorbar(fraction=1, anchor=(0,0), extend='both')
    #~ ax = plt.axes([0.35, 0.25, 0.8, 1], visible=False)
    #~ plt.colorbar(ax=ax, orientation='horizontal', anchor = (0.2, 0), shrink=0.5, pad=0.2, fraction=0.30, drawedges=False, aspect=10)
    plt.figtext(0.01 , 0.01 ,os.path.abspath(pathname) + '/' + inspect.getfile(inspect.currentframe())) # OR os.path.realpath(__file__)
    plt.savefig(FIGDIR+name+'.png')
    plt.close('all')

#~ plot_tll(cflx, 'cflx')

def plot_summer_winter_means(year_stack, name, lat=lat, lon=lon):
	'''
	Plotting function to plot mean summer 
	(JFM) and winter (JAS) climatologies.
	'''
	time_end = np.shape(year_stack)[1]; print time_end #~ plt.figure(figsize=(15, 15), dpi=300)
	cmap = cm.RdBu_r
	curr_map = Basemap(projection='cyl', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='c', area_thresh=100.)
	#~ curr_map = Basemap(projection='spstere',boundinglat=bounding_lat,lon_0=lon_0,resolution='h', area_thresh=50.)
	x, y = curr_map(lon, lat)
	#~ plt.figure(figsize=(12, 8), dpi=50)
	summer = year_stack[:, 0:time_end/4, :, :]
	summer_climatological = np.mean(np.mean(summer, axis=1), axis = 0)
	winter = year_stack[:, (time_end/4)*2:(time_end/4)*3, :, :]
	winter_climatological = np.mean(np.mean(winter, axis=1), axis = 0)
	plt.subplot(211)
	#summer, winter = summer_and_winter_climatological(array, 5)
	#~ vmin = min(np.min(summer), np.min(winter))
	#~ vmax = max(np.max(summer), np.max(winter))
	#~ plt.pcolormesh(summer_climatological, vmin=vmin, vmax=vmax); 
	im = curr_map.pcolormesh(x, y, summer_climatological, cmap=cmap, vmin=vmin, vmax=vmax); plt.axis('tight'); # , vmin=vmin, vmax=vmax
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')
	plt.yticks([-78+16, -78+16+11, -78+16+11+9, -78+16+11+9+7], s2n[20:70:10])
	plt.xticks( np.arange(-160, 180, 60), e2w[20:-1:60])
	#~ plt.colorbar()
	plt.axis('tight')
	#~ plt.title("Simulated summer climatological mean")
	plt.subplot(212)
	#~ plt.pcolormesh(winter_climatological , vmin=vmin, vmax=vmax); 
	im = curr_map.pcolormesh(x, y, winter_climatological,  vmin=vmin, vmax=vmax, cmap=cmap); plt.axis('tight'); # vmin=vmin, vmax=vmax,
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')
	#~ plt.colorbar()
	plt.axis('tight')
	plt.yticks([-78+16, -78+16+11, -78+16+11+9, -78+16+11+9+7], s2n[20:70:10])
	plt.xticks( np.arange(-160, 180, 60), e2w[20:-1:60])
	#~ plt.title("Simulated winter climatological mean")
	ax = plt.axes([0.05, 0.05, 0.9, 0.9], visible=False)
	plt.colorbar(ax=ax, orientation='horizontal', anchor = (0, 0), shrink=0.5, pad=0.25, fraction=0.25, drawedges=False, aspect=55)	
	plt.savefig('/home/nicholas/thesis/data-analysis/output/figures/summer_winter_means_' +  name + '.png')

#~ ys_pco2 = ma.array(np.split(pco2, 10, axis=0))	
plot_summer_winter_means(ma.array(np.split(cflx_new, 10, axis=0)), 'cflx5day')


def plot_stdevs(data, name):
	data /= np.max(np.abs(data), axis=0) 
	year_stack=ma.array(np.split(data, 10, axis=0))
	
	vmin, vmax = 0, 0.5
	#~ plt.figure(figsize=(10, 10), dpi=50)
	curr_map = Basemap(projection='cyl', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='i', area_thresh=100.)
	x, y = curr_map(lon, lat)
	
	plt.subplot(411)
	stdev_all_data = ma.std(data, axis=0)
	im = curr_map.pcolormesh(x, y, stdev_all_data , vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('stdev_all_data'+ longname)
	
	plt.subplot(412)
	annual_means = ma.mean(year_stack, axis = 1)
	stdev_annual_means = ma.std(annual_means, axis=0)
	im = curr_map.pcolormesh(x, y, stdev_annual_means , vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('Standard Deviation of the Annual Averages'+ longname)
	
	plt.subplot(413)
	signal_array = ma.mean(year_stack, axis=0)
	stdev_seasonal = ma.std(signal_array, axis=0)
	im = curr_map.pcolormesh(x, y, stdev_seasonal , vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('stdev_seasonal'+ longname)
	
	plt.subplot(414)
	stdev_all_data = ma.std(data, axis=0)
	signal_array = ma.mean(year_stack, axis=0)
	stdev_seasonal = ma.std(signal_array, axis=0)
	stdev_non_seasonal = stdev_all_data - stdev_seasonal
	#~ stdev_non_seasonal = ma.stdev(noise_array, axis=0)
	im = curr_map.pcolormesh(x, y, stdev_non_seasonal, vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('stdev_non_seasonal' + longname)
		
	plt.savefig('/home/nicholas/masters/figures/newplots/standard_deviations_' + name+ '.png')
	plt.close('all')
	
#~ plot_stdevs(pco2, 'pco2')
#~ plot_stdevs(cflx, 'cflx')

def snr_func(data):
	data /= np.max(np.abs(data), axis=0) 
	year_stack=ma.array(np.split(data, 10, axis=0))
	stdev_all_data = ma.std(data, axis=0)
	signal_array = ma.mean(year_stack, axis=0)
	stdev_seasonal = ma.std(signal_array, axis=0)
	stdev_non_seasonal = stdev_all_data - stdev_seasonal
	return stdev_seasonal/stdev_non_seasonal
	
def snr(cflx, pco2):
	#~ plt.figure(figsize=(10, 10))
	curr_map = Basemap(projection='cyl', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat, resolution='i', area_thresh=100.)
	x, y = curr_map(lon, lat)
	vmin=0
	vmax=10
	V = np.arange(0, 10)
	
	plt.subplot(211)
	im = curr_map.pcolormesh(x, y, snr_func(cflx) , vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('Cflux SNR')
	
	plt.subplot(212)
	im = curr_map.pcolormesh(x, y, snr_func(pco2), vmin=vmin, vmax=vmax, cmap=cmap)
	plt.axis('tight')
	plt.colorbar()
	curr_map.drawcoastlines()
	curr_map.fillcontinents(color='grey',lake_color='aqua')	
	#~ plt.title('Dpco2 SNR')	
	
	plt.savefig(FIGDIR + 'snrs.png')
	plt.close('all')
	
#~ snr(cflx, pco2)

def ts_ft(data, name):
	data /= np.max(np.abs(data), axis=0) 
	plt.subplot(311)
	plt.plot(ma.mean(ma.mean(data, 2), 1))
	plt.axis('tight')
	plt.subplot(312)
	plt.plot(signal.detrend(ma.mean(ma.mean(data, 2), 1), axis=0))
	plt.axis('tight')
	plt.subplot(313)
	plt.plot(np.log(np.fft.rfft(ma.mean(ma.mean(signal.detrend(data, axis=0), 2), 1)))[0:365])
	plt.axis('tight')
	plt.savefig(FIGDIR+'time_freq_'+name+'.png')
	plt.close('all')
	
#~ ts_ft(cflx/30, 'cflx')
#~ ts_ft(pco2, 'pco2')

def ft_2d(data, name):
	#~ data /= np.max(np.abs(data), axis=0) 
	vmin, vmax = 0, 15
	plt.subplot(121)
	ft_lat = np.fft.rfft2(ma.mean(data, 2))[0:np.shape(data)[0]/2, :]
	plt.pcolormesh(np.log(ft_lat), vmin=vmin, vmax=vmax); plt.axis('tight'); plt.colorbar()
	plt.subplot(122)
	ft_lon = np.fft.rfft2(ma.mean(data, 1))[0:np.shape(data)[0]/2, :]
	plt.pcolormesh(np.log(ft_lon), vmin=vmin, vmax=vmax); plt.axis('tight'); plt.colorbar()
	#~ ax = plt.axes([0, 0, 0, 0], visible=False)
	#~ plt.subplot((223), visible=False)
	#~ plt.colorbar(orientation='horizontal', fraction=0.5)	
	plt.savefig(FIGDIR+'ft_'+name+'.png')
	plt.close('all')
	
#~ ft_2d(cflx, 'cflx')
#~ ft_2d(pco2, 'pco2')

def get_boundaries(file_name):
    fname = open(file_name)
    boundary = np.genfromtxt(fname)
    #x = boundary[:, 0]
    #y = boundary[:, 1]
    return boundary
    
'''
ma.mean(cflx, 0)
#plot annual means
ma.mean(ys_pco2, 1)
ma.mean(ys_cflx, 1)
#plot summer winter w takahashi
s_cflx
w_cflx
s_pco2
w_pco2
#plot stdevs
ma.std(pco2, 0)
ma.std(ma.mean(ys_pco2), 1)
ma.std(pco2_sig, 0)
ma.std(pco2_noi, 0)
#plot snrs cflx dpco2
snr_pco2
snr_cflx
#plot timeseries, detrended and fourier
ma.mean(ma.mean(pco2, 2), 1)
ma.mean(ma.mean(det_pco2, 2), 1)
np.fft.rfft(ma.mean(ma.mean(det_pco2, 2), 1)) # log?
#plot 2dfts 1st quad
np.fft.rfft2(ma.mean(pco2, 1))
np.fft.rfft2(ma.mean(pco2, 2))

#Then the brute-force method
?exhaustive search?
'''

'''
def more_stuff(input_array):
    year_stack = np.split(input_array, 10, axis=0)
    normalized_array /= np.max(np.abs(input_array, axis=0)
    detrended_array = detrend(input_array, axis=0)
    annual_cycle = ma.mean(year_stack, axis=0)
    signal_array = np.tile((annual_cycle), (10, 1, 1))
    noise_array = input_array - signal_array
    summer = year_stack[:, 0:90, :, :]
    summer_climatological = np.mean(np.mean(summer, axis=1), axis = 0)
    winter = year_stack[:, 180:270, :, :]
    winter_climatological = np.mean(np.mean(winter, axis=1), axis = 0)
    stdev_signal_array = ma.std(signal_array, 0)
    stdev_noise_array = ma.std(noise_array, 0)
    stdev_all_data = ma.std(input_array, 0)
    snr_lenton_time_domain = stdev_signal_array/stdev_noise_array
    time_lat_array = ma.mean(input_array, 2)
    time_lon_array = ma.mean(input_array, 1)
    ft_time_lat_array = []
    ft_time_lon_array = []
'''





