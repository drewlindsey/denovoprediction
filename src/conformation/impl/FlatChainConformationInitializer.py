class FlatChainConformationInitializer(AbstractConformationInitializer):
	"""Initializes a conformation with linear, 180 degree, angles
	
	Attributes:
		sequence: The sequence to initialize
	"""
	def __init__(self, sequence):
		"""Inits the conformation initializer."""
		self.sequence = sequence
	
	def generateConformation(self):
		"""Generates a straight-line backbone conformation."""
		raise Exception('Not implemented yet, cuz lazy')