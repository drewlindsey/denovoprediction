# TODO move conformation initialization to this class ?
class Conformation(object):
	"""A representation of an amino acid tertiary structure.

	Attributes:
		sequence: the amino acid sequence characters
		conformationInitializer: An object of type AbstractConformationInitializer to create the initial "backbone" (list of Residues)
		conformation: A list of Residue values	
	"""
	def __init__(self, sequence, conformationInitializer):
		"""Inits Conformation with the conformationInitializer"""
		self.sequence = sequence
		self.conformationInitializer = conformationInitializer
		self.conformation = self.conformationInitializer.generateConformation()
	
	def reset(self):
		"""Resets back to the initial conformation using the original ConformationInitializer"""
		self.conformation = self.conformationInitializer.generateConformation()
		
	def reset(self, conformationInitializer):
		"""Resets to an initial conformation using the given ConformationInitializer"""
		self.conformationInitializer = conformationInitializer
		self.conformation = self.conformationInitializer.generateConformation()
		
	def get(self, position):
		"""Returns a Residue at the given position"""
		return self.conformation[position]