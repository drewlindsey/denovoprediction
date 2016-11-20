class Fragment(object):
	def __init__(self, residues, k):
		self.k = k
		self.residues = residues
	
	def getResidues(self):
		return self.residues
		
	def getResidue(self, index):
		return self.residues[index]
		
	def getK(self):
		return self.k