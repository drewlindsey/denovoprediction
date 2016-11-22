class AbstractConformationInitializer(object):
	"""An abstract base class for conformation initialization. Determines
	how the conformation is initialized, e.g. straight-line, 'nearly' straight-line,
	or a 'match' from the PDB
	
	Attributes:
		sequence: The sequence to initialize into a list of Residue objects
	"""
	def __init__(self, sequence):
		"""Inits the Conformation Initializer"""
		raise Exception('Abstract class has no implementation')
	
	def generateConformation(self):
		"""Generates the initialized conformation using this class' methods"""
		raise Exception('Abstract class has no implementation')