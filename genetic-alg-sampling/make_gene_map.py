# ~/Desktop/ga-session-2012-11-12/make_gene_map.py

'''Make a gene map that maps a binary string
to a coordinate of an array as well as a value 
for that coordinate.'''

import numpy as np
import itertools
import load_cflux_masked
from numpy import ma

class GeneMap(object):

    def __init__(self):
        self.gene_map = {}
        self.location_dict = {}
        self.actual_data_dict = {}
        self.location_dict_stdevs = {}
        self.array = []
        self.time_len = 0
        self.lat_len = 0
        self.lon_len = 0
        self.last_valid_binary_string = ''
        self.string_length = 0
        self.count = 0
        self.what_data = ''

    def get_array_attributes(self):
        lat_end = 45
        lon_end = 180
        ### 5 day data or daily data
        time_end = 730
        #time_end = 3650
        ### Choose array (decadel mean, annual mean or all data)
        ### All Data
        self.what_data = 'AllData'
        self.array = load_cflux_masked.load_file(time_end = time_end, lat_end=lat_end, lon_end=lon_end)
        ### decadel mean
        #~ self.what_data = 'DecadalMean'
        #~ self.array = ma.mean(self.array, axis=0)
        #~ self.array = np.reshape(self.array, (1, lat_end, lon_end))
        ### annual mean
        self.what_data = 'AnnualCycle'
        self.array =  ma.mean(np.split(self.array, 10, axis=0), axis=0)
        ###
        self.array_shape = np.shape(self.array)
        print self.array_shape
        self.count_non_masked = ma.count(self.array)
        self.time_len = self.array_shape[0] # need to set interpolated and masked array time_end to be equal NB!!!
        self.lat_len = self.array_shape[1]
        self.lon_len = self.array_shape[2]
        self.string_length = len(bin(self.count_non_masked)[2:])
        for item in itertools.product(range(self.lat_len), range(self.lon_len)):
            self.actual_data_dict[item] = np.std(self.array[: , item[0], item[1]])


    def make_gene_map_2(self):
        '''
        The method that takes the attributes from the array and uses
        them to create a gene map for the array.
        The gene map is a dictionary which has a binary string as a key.
        The binary string is created by creating a binary bit string of 
        an appropriate length.
        The length is calculated 
        '''
        count = 0
        self.iterator_one = itertools.product(range(self.time_len), range(self.lat_len), range(self.lon_len))
        ### Assign a binary string a location and a value from the data
        print "\n"
        print "Creating gene-map dictionary... \n"
        print "Assigning valid locations to binary strings! \n"
        for x_valid in self.iterator_one:
            binary_string = bin(count)[2:]
            while len(binary_string) < self.string_length: # removed minus one (-1) NB
                binary_string = '0' + binary_string
            self.gene_map[binary_string] = {}
            if ma.is_masked(self.array[x_valid]):
                pass
            else:
                self.gene_map[binary_string]['coordinate'] = tuple(x_valid)
                self.gene_map[binary_string]['value'] = self.array[x_valid]
                self.location_dict[x_valid[1:3]] = []
                self.location_dict_stdevs[x_valid[1:3]] = 0
                count += 1
        self.last_valid_binary_string = binary_string
        binary_string_old = binary_string
        not_valid_first = int(binary_string, 2) + 1
        not_valid_last = int('1'*(self.string_length), 2) # added minus one just for nonmasked version NB
        self.count = count
        
        if self.count == self.count_non_masked:
            print "The counter corresponds with the non-masked count! \n"
        ### Pad the dictionary to give binary strings some value
        print "Assigning left over binary strings to non-existant locations! \n"
        self.iterator_two = itertools.product(range(self.time_len), range(self.lat_len), range(self.lon_len))
        #~ for x_not_valid in range(not_valid_first, not_valid_last+1):
        count_2 = not_valid_first
        for x_not_valid in self.iterator_two:
            #~ binary_string = bin(x_not_valid)[2:] # DOES IT NEED TO BE PADDED
            binary_string = bin(count_2)[2:]
            while len(binary_string) < self.string_length: # removed minus one (-1) NB
                binary_string = '0' + binary_string
            self.gene_map[binary_string] = {}
            self.gene_map[binary_string]['coordinate'] = (999, 999, 999) #x_not_valid
            self.gene_map[binary_string]['value'] = 1E09 #self.array[x_valid]
            if count_2 == not_valid_last:
                break
            else:
                count_2 += 1
        print "There are %d valid locations. \n"  %count
        print "The last binary string is: ", binary_string
        print "The last binary string assigned to a valid locations is :", binary_string_old
        print "The length of binary string is: ", self.string_length
        print "The non-valid locations fall between %d and %d. \n" %( not_valid_first, not_valid_last)
        print "Is the array masked?: \n", ma.isMA(self.array)
        print "The gene-map has been created! \n"


			
