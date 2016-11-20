class AbstractConformationSampler(object):
	def __init__(self, initialConformation, seefModel, fragLib):
		self.conformation = initialConformation
		self.seef = seefModel
		self.fragLib = fragLib

	def next(self):
		raise Exception('Abstract class has no implementation')
		
	def hasNext(self):
		raise Exception('Abstract class has no implementation')
		
	def minimum(self):
		raise Exception('Abstract class has no implementation')