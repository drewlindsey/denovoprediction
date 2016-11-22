class AbstractPipeline(object):
	"""An abstract base class for a de novo pipeline. An implementation should have all 
	aspects of the pipeline to generate the predicted structure.
	
	TODO: pass in the required objects to build this pipeline??
	
	Attributes:
		sequence: The amino acid sequence
	"""
	def __init__(self, sequence):
		"""Inits the pipeline with the given sequence"""
		raise Exception('Abstract class has no implementation')
		
	def generateDeNovoPrediction():
		"""Execute the pipeline and obtain the predicted structure conformation"""
		raise Exception('Abstract class has no implementation')