class Conformation(object):
	def __init__(self, sequence, backboneInitializer):
		self.sequence = sequence
		self.backboneInitializer = backboneInitializer
		
	def initialize():
		self.conformation = self.backboneInitializer.initialize(self.sequence)
		
	def reset(self, backboneInitializer):
		self.backboneInitializer = backboneInitializer
		self.conformation = self.backboneInitializer.initialize(self.sequence)
		
	def get(position):
		return self.conformation[position]