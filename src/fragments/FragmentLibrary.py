class FragmentLibrary(object):
	def __init__(self, sequence, frag3, frag9):
		self.sequence = sequence
		self.fragments3 = frag3
		self.fragments9 = frag9
		
	def get3merFragments(self):
		return self.fragments3
		
	def get9merFragments(self):
		return self.fragments9
		
	def get3merFragment(self, position):
		return self.framents3[position]
		
	def get9merFragment(self, position):
		return self.fragments9[position]