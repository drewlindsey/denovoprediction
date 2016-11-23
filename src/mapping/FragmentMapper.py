class FragmentMapper(object):
	"""A ROBETTA -> Fragment mapper. Takes in a ROBETTA text file and parses it creating all the 
	fragment lists.
	
	Attributes:
		k: The k-mer length (k == 3 or k == 9)
	"""
	def __init__(self, k):
		"""Inits the FragmentMapper"""
		self.k = k
		
	def mapFileToFragments(self, inputFileName):
		"""Takes the input file and outputs the fragment list"""
		fragments = []
		start = False
		count = 0
		currFragCount = 0
		with open(inputFileName) as robFile:
			for line in robFile:
				if count == 2 and start:
					start = True
					count = 0
				if count == 4 and not start:
					count = 0
					if currFragCount == 200:
						currFragCount = 0
				if count < 4 and not start:	
					currFragCount += 1
					data = line[18:44]
					print data
				count += 1
				
		return fragments