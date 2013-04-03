# rbf_using_fortran.py

import numpy as np
from FModule import *
from scipy import sparse,linalg
from scipy.sparse import linalg

# TODO: Question: Why are there so many zeros in the 'interpolated' data?
# values = values at each location for each year
# locations = the locations (same for each year)
# valid_locations (There are about 591300 each year)
# data = the data that is used to choose locations (year average)
# s_out has what shape????

# paste/run Main.py
# paste/run process_locations.py
# then run rbf_using_fortran.py
# then call the functions with correct inputs

# formatted_list = format_lists(lat_lon_loc, values)
# out_list = rbf_on_each_year(formatted_list[1], formatted_list[0], locations)
# interpolated_data = interp_each_year(data, out_list)

# make lat_lon_loc and values into same format as can be used in FModule
# ie. lat_lon_loc is array with shape (6903, 3) and values is list with shape (69030)
def format_lists(sample_locations, values):
    _sample_locations = [0.0, 0.0, 0.0]
    for item in sample_locations:
        for el in sample_locations[item]:
            _sample_locations = np.vstack((_sample_locations, np.append(item, el)))
    _sample_locations = _sample_locations[1:, :]
    _values = []
    for item in values:
        for el in values[item]:
            _values.append(el)
    _values_split = np.split(np.array(_values), 10)
    return _sample_locations, _values_split
    
# for each year do a RBF using FModule and return the 
# interpolated array

def rbf_on_each_year(sample_values, sample_locations, locations):
    r = 20.0
    out_list = []
    for year in range(10):
        print year
        r_vector,mat = fortranmodule.rbf_function_matrix(r,sample_values[year],sample_locations,locations)
        r1 = sparse.linalg.spsolve(mat.copy(),r_vector.copy())
        s_out = fortranmodule.rbf_evaluate(r,sample_values[year],sample_locations,locations,r1)[:,0]
        out_list.append(s_out)
    return out_list
    
def data_interpolater(data, s_out):
    data_interp = np.zeros((data.shape[0],data.shape[1],data.shape[2]))
    for i in range(s_out.shape[0]):
        t = locations[valid_locations[i],0]
        lon = locations[valid_locations[i],1]
        lat = locations[valid_locations[i],2]
        data_interp[t,lon,lat] = s_out[i]
    return data_interp

def interp_each_year(data, out_list):
    interpolated_data = np.zeros((10, 73, 45, 180))
    for year in range(10):
        print year
        interpolated_data[year, :, :, :] = data_interpolater(data, out_list[year])
    return interpolated_data


def calculate_errors():
    for year in range(10):
        error[year, :, :, :] = data[year, :, :, :] - interpolated_data[year, :, :, :]

# The End       
