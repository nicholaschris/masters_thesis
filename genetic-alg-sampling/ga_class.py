#! usr/bin/python
# ~/Desktop/ga-session-2012-11-12/ga_class.py

'''All the attributes and methods for the ga class.'''

import numpy as np
import random
from numpy import ma
from operator import itemgetter
from scipy.interpolate import Rbf
import datetime
import make_gene_map
import make_population
import os
import itertools

class GeneticAlgorithm(make_gene_map.GeneMap, make_population.InitPopulation):

    def __init__(self):
        make_gene_map.GeneMap.__init__(self)
        make_population.InitPopulation.__init__(self)
        self.sum_fitness = 0
        self.new_pop = []
        self.least_fitness = 0
        self.fittest_fitness = 0
        self.mutation_rate = 0.000001
        self.diff_array = []
        self.mask = []
        self.xxx, self.yyy, self.zzz = 0, 0, 0
        self.array_mean = 0
        self.array_stdev = 0
        self.array_range = 0

    def getSolutionCosts(navigationCode):
        '''
        This should go in fitness function.
        '''
        fuelStopCost = 15
        extraComputationCosts = 8
        thisAlgorithmBecomingSkynetCost = 999999999
        waterCrossingCost = 45
     
    def get_mask(self):
        self.array_mean = ma.mean(self.array)
        self.array_stdev = ma.std(self.array)
        self.array_range = ma.max(self.array) - ma.min(self.array)
        print "The mean is %f, the stdev is %f, the range is %f." %(self.array_mean, self.array_stdev, self.array_range)
        from scipy.io.netcdf import netcdf_file as NetCDFFile
        ### get landmask
        nc = NetCDFFile(os.getcwd()+ '/../data/netcdf_files/ORCA2_landmask.nc','r')
        self.mask = ma.masked_values(nc.variables['MASK'][:, :self.time_len, :self.lat_len, :180], -9.99999979e+33)
        nc.close()
        self.xxx, self.yyy, self.zzz = np.lib.index_tricks.mgrid[0:self.time_len, 0:self.lat_len, 0:180]

    def calc_chrom_rbf_part_1(self, index):
        '''
        Calculates the fitness of a chromosome according to 
        the square of the difference of an RBF interpolation of
        the sampled points from the chromosome and the
        model data.
        Used in add_fitness_key method.
        '''
        chromosome = self.population[index]['chrom_list']
        #~ print "Length of chromosome: ", len(chromosome) - self.string_length
        gene_stepper = 0
        values_list = []
        coord_list = []
        gene_length = self.string_length # hard code ? AGAIN -1 ???
        while gene_stepper < len(chromosome) - gene_length: # The gene stepper goes from gene to gene o chromosome
            gene_list = chromosome[gene_stepper:gene_stepper+self.string_length]
            current_gene = '' # initialises current gene
            bit_stepper = 0 # counter for the bits
            while bit_stepper < gene_length:
                current_gene = current_gene + gene_list[bit_stepper] #constructs the current gene
                bit_stepper += 1
            values_list = np.append(values_list, self.gene_map[current_gene]['value']) #appends the value of the current gene to a list
            coord_tuple = self.gene_map[current_gene]['coordinate'] # assigns  current coordinate to a variable
            coord_list.append(coord_tuple)
            gene_stepper += gene_length # goes to the next gene in the loop
        coords_vals = dict(zip(coord_list, values_list))
        no_of_locations = np.size(coords_vals.keys())/3 # Divide by 3 really???? Yes np.size != len
        #print "Coords_vals.keys() has length: ", no_of_locations*3
        if coords_vals.has_key((999, 999, 999)):
            #~ print "Invalid location selected by individual!"
        #~ print coords_vals
            self.invalid_count+=1
        return values_list, coords_vals, no_of_locations
        #return coords_vals

    def calc_chrom_rbf_part_2(self, index, coords_vals):
        x_nodes = []	
        y_nodes = []
        z_nodes = []
        values_list = []
        for item in coords_vals[1]:
            #print item
            x_nodes.append(item[0])
            y_nodes.append(item[1])
            z_nodes.append(item[2])
            values_list = np.append(values_list, coords_vals[1][item]) # What is the difference between the old and new values list
        ### Trying out 3d RBF (insert into ga_class)
        ### interpolate sample data - use RBF
        #~ rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)
        rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)
        ZI = rbf(self.xxx, self.yyy, self.zzz)
        ZI = ZI*self.mask[0, :, :, :]
        self.array = self.array*self.mask[0, :, :, :]
        GI = self.array[:, :, :]
        #~ self.diff_array = (GI - ZI)**2
        fitness = np.sqrt(np.mean((GI - ZI)**2))
        self.population[index]['fitness'] = fitness
        
    def calc_chrom_fast(self, index, coords_vals):
        self.population[index]['fitness'] = \
        np.abs(self.array_mean - ma.mean(coords_vals[0])) + \
        np.abs(self.array_stdev - ma.std(coords_vals[0])) + \
        np.abs(self.array_range - (ma.max(coords_vals[0])-ma.min(coords_vals[0])))/10  + \
        np.abs((self.chromosome_size-1) - coords_vals[2]) #locations
        #~ print "Chromosome size: ",self.chromosome_size
        print "Number of locations is: ", coords_vals[2]
        #~ print "The sample range is: %g. The array range is: %g " % ((ma.max(coords_vals[0])-ma.min(coords_vals[0])), self.array_range)
        #~ print np.abs(self.array_mean - ma.mean(coords_vals[0])), np.abs(self.array_stdev - ma.std(coords_vals[0])), np.abs(self.array_range - (ma.max(coords_vals[0])-ma.min(coords_vals[0])))
        #~ print ma.mean(coords_vals[0]), ma.std(coords_vals[0]), (ma.max(coords_vals[0])-ma.min(coords_vals[0]))
        #~ print "Fitness is: ", self.population[index]['fitness']

    def calc_chrom_stdev(self, index, coords_vals):
        e_squared = []
        for item in self.location_dict:
            self.location_dict[item] = []
        for item in self.location_dict_stdevs:
           self.location_dict_stdevs[item] = 0
        for item in coords_vals[1]:
            if item != (999, 999, 999):
                self.location_dict[item[1:3]].append(coords_vals[1][item]) #
            else:
                pass
        for item in itertools.product(range(self.lat_len), range(self.lon_len)):
            try:
                self.location_dict_stdevs[item] = np.std(self.location_dict[item])
                e_squared.append((self.location_dict_stdevs[item] - self.actual_data_dict[item])**2)
            except KeyError:
                pass

        rmse = np.sqrt(np.mean(e_squared))
        self.population[index]['fitness'] = rmse
        
    def add_fitness_key(self):
        '''
        Goes through each chromosome in the population and adds 
        a fitness key and a value for that key.
        Used in sort_pop method.
        '''
        self.invalid_count=0
        for i in range(len(self.population)):
            coords_vals = self.calc_chrom_rbf_part_1(i)
            self.calc_chrom_rbf_part_2(i, coords_vals)
            #~ self.calc_chrom_fast(i, coords_vals)
        #~ print datetime.datetime.now()
        print "Invalid count is: ",self.invalid_count
        
    def add_fitness_fast(self):
        '''
        Goes through each chromosome in the population and adds 
        a fitness key and a value for that key.
        Used in sort_pop method.
        '''
        self.invalid_count=0
        for i in range(len(self.population)):
            coords_vals = self.calc_chrom_rbf_part_1(i)
            self.calc_chrom_fast(i, coords_vals)
        print "Invalid count is: ",self.invalid_count

    def add_fitness_rmse_stdev(self):
        '''
        Goes through each chromosome in the population and adds 
        a fitness key and a value for that key.
        Used in sort_pop method.
        '''
        self.invalid_count=0
        for i in range(len(self.population)):
            coords_vals = self.calc_chrom_rbf_part_1(i)
            self.calc_chrom_stdev(i, coords_vals)
        print "Invalid count is: ",self.invalid_count

    def sort_pop(self):
        '''
        Uses the fitness key of the chromosomes
        in the poputlation to sort the population form 
        fit to least fit.
        '''
        self.population = sorted(self.population, key=itemgetter('fitness')) 
        self.fittest_fitness = self.population[0]['fitness']
        self.least_fitness = self.population[-1]['fitness']
        print "The fittest individual has a fitness of %g. The least fit individual has a fitness of %g" % (self.fittest_fitness, self.least_fitness)
        print '*'*79

    def tourn_sel(self):
        '''
        Chooses two random individuals from the
        population.  
        The individual with the lowest value for the fitness key
        becomes a parent for the next generation.
        '''
        x = random.randint(0, 19)
        player1 = self.population[x]
        y = random.randint(0, 19)
        player2 = self.population[y]
        if player1['fitness'] <= player2['fitness']:
            parent = player1['chrom_list']
        else:
            parent = player2['chrom_list']
        return parent
        #~ print x, y
    
    def select_parent_from_tournament(self):
        #~ '''
        #~ Selects a parent using tournament selection.
        #~ '''
        #~ return self.tourn_sel()
        '''
        Chooses two random individuals from the
        population.  
        The individual with the lowest value for the fitness key
        becomes a parent for the next generation.
        '''
        x = random.randint(0, 19)
        player1 = self.population[x]
        y = random.randint(0, 19)
        player2 = self.population[y]
        if player1['fitness'] <= player2['fitness']:
            parent = player1['chrom_list']
        else:
            parent = player2['chrom_list']
        return parent

    def crossover(self):
        '''
        Selects two parents from the population
        and then mutates the parents
        and then creates two children from 
        the two parents using crossover.
        '''
        #~ crossover_point = random.randint(0, len(self.population[0]['chrom_list']))
        crossover_point = random.randint(0, self.chromosome_size)*(self.string_length)
        #~ print len(self.population[0]['chrom_list'])
        #~ print self.string_length*200
        #~ assert crossover_point%self.string_length == 0
        #~ print "Crossover on gene: ", crossover_point
        parent1 = self.select_parent_from_tournament()
        parent2 = self.select_parent_from_tournament()
        parent1 = self.mutate(parent1)
        parent2 = self.mutate(parent2)
        child1 = parent1[:crossover_point] + parent2[crossover_point:]
        child2 = parent1[crossover_point:] + parent2[:crossover_point]
        return child1, child2

    def mutate(self, chromosome):
        chromosome = list(chromosome)
        for i in range(len(chromosome)):
            if random.random() < self.mutation_rate:
                print 'mutation on %i' % i
                print chromosome[i]
                if chromosome[i] =='0':
                    chromosome[i] = '1'
                else:
                    chromosome[i] = '0'
        return chromosome

    def make_new_pop(self):
        '''
        Uses the crossover function ten times
        in order to get twenty children
        and make a new population.
        '''
        self.new_pop = []
        for i in range(10):
            dictionary1= {}
            dictionary2 = {}
            dictionary1['chrom_list'], dictionary2['chrom_list'] = \
            self.crossover()
            #dictionary1['fitness'], dictionary2['fitness'] = 0, 0
            self.new_pop = np.append(self.new_pop, [dictionary1, dictionary2])

    def elitism(self):
        '''
        Selects a random individual from the new population
        and replaces it with the fittest.
        '''
        r = random.randint(0, 19)
        #print self.population[0]
        self.new_pop[r] = self.population[0]
        ### Double elitism
        #r = random.randint(0, 19)
        #self.new_pop[19] = self.population[1]


