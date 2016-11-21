class AbstractConformationInitializer(object):
	def __init__(self, sequence):
		raise Exception('Abstract class has no implementation')
	
	def generateConformation(self):
		raise Exception('Abstract class has no implementation')