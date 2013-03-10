#!usr/bin/python
# ~/Desktop/ga-session-2012-11-12/ga.py
# modified 2012-11-12-19:27

import datetime
import numpy as np
import make_population
import make_gene_map
import matplotlib
import ga_class
import matplotlib.pyplot as plt
from numpy import ma
from scipy.interpolate import Rbf
from matplotlib import cm
import cPickle
import os
from scipy.io.netcdf import netcdf_file

output_dir = os.getcwd() + '/output/plots/'

### Some stuff to save time and date on plot and filenames
dt = str(datetime.datetime.now())
date = dt[:10]
time = dt[11:16]

genetic = ga_class.GeneticAlgorithm() 		
genetic.get_array_attributes()
genetic.make_gene_map_2()
genetic.get_gene_size(genetic.string_length)
genetic.set_up_new_pop()
print "The time axis has length: ",genetic.time_len
print "The chromosome size is: ", genetic.chromosome_size
genetic.get_mask()

fittest_list = []		# Use the fittest list to plot evolution of algorithm DOES IT WORK?
least_fit_list = []	# likewise for least fit list

print datetime.datetime.now()
user_input = raw_input("The GA is set up, press y to continue.")
if user_input == 'y':
	print "Continuing"
	
# Set number of generations here:
NUM_GEN = 1000
print "Running Genetic Algorithm for %d generations." % NUM_GEN
gen_count = 0
while gen_count < NUM_GEN:
    #~ genetic.add_fitness_key()
    genetic.add_fitness_fast()
    #~ genetic.add_fitness_rmse_stdev()
    fittest_list = np.append(fittest_list, genetic.fittest_fitness)
    least_fit_list = np.append(least_fit_list, genetic.least_fitness)
    genetic.sort_pop()
    genetic.make_new_pop()
    genetic.elitism()
    #~ print genetic.population[0]['fitness']
    #~ print genetic.new_pop[0]['fitness']
    genetic.population = genetic.new_pop
    print "There are %g generations left" %(NUM_GEN-gen_count) 
    gen_count+=1

print datetime.datetime.now()
plt.close('all')

### A plot to visualise how the GA evolved ###
plt.ylim(ymin=0, ymax = 400)
plt.plot(fittest_list, label='Fittest Individuals')
plt.plot(least_fit_list, label='Least Fit Individuals')
plt.title('Genetic Algorithm - ' + date + time)
plt.xlabel('Generation')
plt.ylabel('fitness')
plt.legend(loc=7)
plt.savefig(output_dir + 'GA' + date + '_' + time + '.png')
### End plotting ###
 
genetic.add_fitness_fast()
genetic.sort_pop() 
solution = genetic.population[0]['chrom_list']

the_fittest_fitness = genetic.population[0]['fitness']

### This saves the coordinate tuples in a list ###
def calc_chrom(chromosome, dictionary):
    '''
    This saves the coordinate tuples in a list.
    '''
    gene_stepper = 0
    values_list = []
    coord_list = [(0, 0, 1)]
    gene_length = genetic.string_length
    while gene_stepper < len(chromosome) - gene_length:
        gene_list = chromosome[gene_stepper:gene_stepper+genetic.string_length]
        current_gene = ''
        bit_stepper = 0
        while bit_stepper < gene_length:
            current_gene = current_gene + gene_list[bit_stepper]
            bit_stepper += 1
        values_list = np.append(values_list, dictionary[current_gene]['value'])
        coord_tuple = dictionary[current_gene]['coordinate']
	coord_list.append(coord_tuple)
	gene_stepper += gene_length
    coord_list.pop(0)
    print "Size of Coord_list: ", np.size(coord_list) # DEBUG
    return coord_list, values_list



coord_list, values_list = calc_chrom(solution, genetic.gene_map)
coords_vals = dict(zip(coord_list, values_list))
no_of_locations = np.size(coords_vals.keys())
print "Coords_vals.keys() has length: ", no_of_locations
if coords_vals.has_key((999, 999, 999)):
    print "Popping invalid location!"	
    coords_vals.pop((999,999,999))
    


# Print the list of coordinates to a file!
pop_size = str(np.size(values_list))
data_used = genetic.what_data
fitness_function_used = 'stats' # 'rbf' or 'stats'
filename = os.getcwd() +'/output/coords/'+'coords-'+str(the_fittest_fitness)+data_used+'-'+date+'_'+time+\
'-'+str(pop_size)+'-'+str(NUM_GEN)+'-'+fitness_function_used+'.pkl'
f= open(filename, 'wb')
cPickle.dump(coords_vals, f)
f.close()

print 'Coordinates saved to ' + filename
print genetic.array_shape

'''
print coord_list
print genetic.array_mean
print ma.mean(values_list)

x_nodes = []
y_nodes = []
z_nodes = []
values_list = []
for item in coords_vals:
    print item
    x_nodes.append(item[0])
    y_nodes.append(item[1])
    z_nodes.append(item[2])
    values_list.append(coords_vals[item])

# There is a problem here!!!
print "The number of locations is: ", np.size(values_list)
### get landmask
nc = netcdf_file(os.getcwd() +'/../data/netcdf_files/ORCA2_landmask.nc','r')
mask = ma.masked_values(nc.variables['MASK'][:, :genetic.time_len, :genetic.lat_len, :genetic.lon_len], -9.99999979e+33)
nc.close()

### Trying out 3d RBF (insert into ga_class)
### meshgrid 3d
xxx, yyy, zzz = np.lib.index_tricks.mgrid[0:genetic.time_len, 0:genetic.lat_len, 0:genetic.lon_len]

rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)
interpolated_data = rbf(xxx, yyy, zzz)
interpolated_data = interpolated_data*mask[0, 0:genetic.time_len, :, :]
data = genetic.array[:, :, :]
diff_array = np.sqrt((data - interpolated_data)**2)

plt.close('all')
plt.figure(figsize=(10, 15))
### Plot Interpolated data, difference and simluated data
plt.subplot(3, 1, 1)
plt.pcolormesh(ma.mean(interpolated_data, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
#plt.scatter(z_nodes, y_nodes, len(x_nodes), values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.title('RBF interpolation'); plt.axis('tight'); plt.colorbar()
#~ plt.savefig(output_dir + 'rbf2d' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.subplot(3, 1, 2)
plt.pcolormesh(ma.mean(diff_array, axis=0), cmap=cm.jet, vmin=0, vmax=10); 
#plt.scatter(z_nodes, y_nodes, len(x_nodes), values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.title('Residuals'); plt.axis('tight'); plt.colorbar()
#~ plt.savefig(output_dir + 'diff_array' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
#~ plt.close('all')
plt.subplot(3, 1, 3)
plt.pcolormesh(ma.mean(data, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
plt.scatter(z_nodes, y_nodes, len(x_nodes), values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.title('Simulated data and optimal sampling locations'); plt.axis('tight'); plt.colorbar()
plt.savefig(output_dir + 'ga/ga_output' + '-' + date + '_' + time + '.png')
plt.close('all')

pop_size = str(np.size(values_list))
data_used = genetic.what_data
fitness_function_used = 'stats' # 'rbf' or 'stats'

fitness = np.sqrt(np.mean((data - interpolated_data)**2))
print np.shape(data)
print np.shape(interpolated_data)
print np.mean(data), np.mean(interpolated_data), np.std(data), np.std(interpolated_data)
print "The fitness of the solution is: %g" % fitness # This answer should be the same as the final fittest solution???

'''

