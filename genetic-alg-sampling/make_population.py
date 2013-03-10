# ~/Desktop/ga-session-2012-11-12/make_population.py

'''
A tool to make a population of chromosomes
that correspond to the gene_map created by the
make_gene_map module.

'''


import random
import numpy as np

# Change chromosome_size to change number of locations selected for sampling

class InitPopulation(object):

	def __init__(self, pop_size=20, chromosome_size=10000):
		'''
		Lenton uses approx 640 measurements per year?
		4 x 12 x 13 per year.
		To sample one time step or decadel average:
		12 x 13 = approx 150.
		Sample all the data would use approx 5000?
		'''
		self.population = [] # was {}
		self.pop_size = pop_size
		self.chromosome_size = chromosome_size # is number of locations
		self.gene_size = 0

 
	def get_gene_size(self, number): # possibly not needed
		'''
		In order to make a chromosome, the number
		of bits in each gene of the chromosome needs to be known.
		This is number can be obtained form the genemap.
		'''
		self.gene_size = number 
		# print self.gene_size # debug

	def make_chromosome(self):
		'''
		Creates a chromosome of randomly chosen 0s and 1s
		using the size of the gene * the size of the chromosome.
		The size of the chromosome is the number of locations 
		that can be chosen to do the RBF interpolation.
		
		'''
		self.new_chromosome = []
		for i in range(self.gene_size*self.chromosome_size):
			self.new_chromosome.append(str(random.randint(0,1)))
		# print self.new_chromosome # debug
		return self.new_chromosome

	def make_chromosome_erok81(self):
		new_chromosome = random.getrandbits(self.gene_size * self.chromosome_size)
		return list('{:b}'.format(new_chromosome))
		
	def set_up_new_pop(self):
		'''
		According to the user given population size, a list of chromosomes is created.
		The fitness of each chromosome can be calculated and the list can be sorted.
		'''
		list = [self.make_chromosome() for i in range(self.pop_size)]
		for i in range(self.pop_size):
			dictionary = {}
			dictionary['individual_no'] = i
			dictionary['chrom_list'] = list[i]
			self.population = np.append(self.population, dictionary)
			
	def set_up_new_pop_erok81(self):
		for i in range(self.pop_size):
		    dictionary = {'individual_no':i, 'chrom_list':self.make_chromosome_erok81()}
		    self.population = np.append(self.population, dictionary)
		
		
