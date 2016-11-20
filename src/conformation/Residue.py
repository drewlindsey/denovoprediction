class Residue(object):
	def __init__(self, type, angles):
		self.aminoAcidType = type
		self.angles = angles
	
	def getType(self):
		return self.aminoAcidType
		
	def getAngles(self):
		return self.angles