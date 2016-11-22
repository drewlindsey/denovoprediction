class Fragment(object):
	"""A representation of a k-mer fragment
	
	Attributes:
		residues: The residue list making up this fragment
		k: The number of residues in the fragment
	"""
	def __init__(self, residues, k):
	"""Inits a Fragment with the given residues and the k count of angles"""
		self.k = k
		self.residues = residues
	
	def getResidues(self):
	"""Gets the Residue list for the faagment"""
		return self.residues
		
	def getResidue(self, index):
	"""Gets the residue at the given index"""
		return self.residues[index]
		
	def getK(self):
	"""Gets the k count of residues"""
		return self.k