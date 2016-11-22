import random

class ConformationSampler(AbstractConformationSampler):
	"""An implementation of AbstractConformationSampler. Uses the SEEF and fragment library to generate a new sample and test
	its energy. End goal is to determine the minimum energy conformation. Works similar to an iterator.
	
	Attributes:
		initialConformation: The initial backbone conformation
		seefModel: The SEEF
		fragLib: The k-mer neighbor fragment library
	
	"""
	def __init__(self, initialConformation, seefModel, fragLib):
		"""Inits this conformation sampler so iterations can begin"""
		self.conformation = initialConformation
		self.seef = seefModel
		self.fragLib = fragLib
		self.minimumConformation
		self.kmax = 5000
		self.k = 0
		self.emax = -1
		self.e = self.seef.computeEnergy(self.conformation)
		self.temp = 1000
		self.maxTemp = 1000
		self.minTemp = 10
		self.ebest = e
	
	def next(self):
		"""Generates the next conformation using the metropolis algorithm"""
		dummy = self.conformation
		startPos = random.randint(0, len(dummy))
		
		prob9 = (self.temp - self.minTemp) / (self.maxTemp - self.minTemp)		
		count = 9 if prob9 > random.random() else 3	
		fragment =  self.fragLib.get9merFragment(startPos) if count == 9 else self.fragLib.get3merFragment(startPos)
		for i in range(startPos,startPos+count):
			dummy[i] = fragment[i-startPos]
		
		energy = self.seef.computeEnergy(dummy)
		
		if P(self.e, energy, self.temp) < random.random():
			self.conformation = dummy
			self.e = energy
		
		if(energy < self.e)
			self.minimumConformation = dummy
			self.ebest = energy
		
		self.k = self.k + 1
		self.temp = self.temp - (self.maxTemp - self.minTemp) / k
		return self.conformation
	
	def hasNext(self):
		"""Checks if the conditions have been fulfilled for this sampling."""
		return self.k < self.kmax and self.e > self.emax
	
	def minimum(self):
		"""The current minimum-energy conformation."""
		return	self.minimumConformation	