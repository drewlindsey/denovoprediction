class AbstractConformationSampler(object):
	"""An abstract base class for the ConformationSampler. 
	
	Controls sampling the conformation space and testing the energy. Acts
	as though it is an iterator and can be continually iterated upon until
	the hasNext conditions have been completed.
	
	Attributes:
		initialConformation: The initial backbone conformation
		seefModel: The SEEF
		fragLib: The FragmentLibrary object to sample from
	"""
	
	def __init__(self, initialConformation, seefModel, fragLib):
		"""Inits the ConformationSampler"""
		self.conformation = initialConformation
		self.seef = seefModel
		self.fragLib = fragLib

	def next(self):
		"""Generate the next conformation"""
		raise Exception('Abstract class has no implementation')
		
	def hasNext(self):
		"""Checks if the sampling is complete"""
		raise Exception('Abstract class has no implementation')
		
	def minimum(self):
		"""The current minima conformation"""
		raise Exception('Abstract class has no implementation')