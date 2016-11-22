class FragmentLibraryGenerator(AbstractFragmentLibraryGenerator):
	"""Takes a ROBETTA file and generates the FragmentLibrary object (9-mer and 3-mer fragment lists)
	
	Attributes:
		sequence: The sequence (TODO might only  need the length)
		mapper: Maps a ROBETTA line to a Residue object
	"""
	def __init__(self, sequence, residueMapper):
		"""Inits the FragmentLibraryGenerator"""
		self.sequence = sequence
		self.mapper = residueMapper
		
	def generateLibrary(self):
		"""Generates the FragmentLibrary for the sequence"""
		raise Exception('Not implemented exception')