class FragmentLibrary(object):
	"""A representation of a library of fragments
	
	This library contains both 3-mers and 9-mers which are lists containing at each index
	of the sequence a 3 or 9 residue fragment for that position.
	
	Attributes:
		length: The length of the target sequence
		frag3: The 3-mer fragment lists (2-D array of Fragment objects) (dimensions are typically length x 200)
		frag9: The 9-mer fragment lists (2-D array of Fragment objects) (dimensions are typically length x 200)
	"""
	def __init__(self, length, frag3, frag9):
		"""Initializes the FragmentLibrary for a sequence of length length with the 2D fragment lists"""
		self.sequence = sequence
		self.fragments3 = frag3
		self.fragments9 = frag9
		
	def get3merFragments(self, index):
		"""Gets the 3-mer Fragment 1D list for the given sequence index"""
		return self.fragments3[index]
		
	def get9merFragments(self):
		"""Gets the 9-mer Fragment 1D list for the given sequence index"""
		return self.fragments9[index]
		
	def get3merFragment(self, index, position):
		"""Gets a 3-mer Fragment for the given sequence index at the given position in the fragment list"""
		return self.framents3[index][position]
		
	def get9merFragment(self, index, position):
		"""Gets a 9-mer Fragment for the given sequence index at the given position in the fragment list"""
		return self.fragments9[index][position]