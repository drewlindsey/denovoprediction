class Residue(object):
	"""A representation of a single amino acid and its torsion angles
	
	Attributes:
		aminoAcidType: The type of amino acid
		angles: The torsion angles dictionary {theta, phi, omega}	
	"""
	
	def __init__(self, type, angles):
		"""Inits the Residue with the given type and angle dictionary"""
		self.aminoAcidType = type
		self.angles = angles
	
	def getType(self):
		"""Returns the type of amino acid"""
		return self.aminoAcidType
		
	def getAngles(self):
		"""Returns the dictionary of angles for the Residue. Valid keys are 'theta', 'phi', 'omega'"""
		return self.angles