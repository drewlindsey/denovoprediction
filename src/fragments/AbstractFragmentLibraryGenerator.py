class AbstractFragmentLibraryGenerator(object):
	"""An abstract base class for generating the FragmentLibrary
	
	Attributes:
		sequence: The sequence used to generate the nearest neighbor fragment library
		residueMapper: The residue mapping class to generate residues from the ROBETTA library (or some other source)
	"""
	
	def __init__(self, sequence, residueMapper):
		"""Inits the fragment library generator"""
		raise Exception('Abstract class has no implementation')
		
	def generateLibrary(self):
		"""Generates the FragmentLibrary for the classes sequence"""
		raise Exception('Abstract class has no implementation')