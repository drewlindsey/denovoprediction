class AbstractSeef(object):
	"""An abstract class template for the Statistical Effective Energy Function (SEEF)
	
	Attributes:
		options: the options for the SEEF function (if any)
	"""
	
	def __init__(self, options=None):
		"""Inits the SEEF with the given (optional) options"""
		raise Exception('Abstract class has no implementation')
		
	def computeEnergy(self, conformation):
		"""Computes the energy of the given conformation as described by the implementation details"""
		raise Exception('Abstract class has no implementation')