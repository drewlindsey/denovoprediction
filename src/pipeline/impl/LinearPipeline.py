class LinearPipeline(AbstractPipeline):
	"""A linear backbone initialized default pipeline.
	
	Attributes:
		sequence: The sequence to determine the tertiary structure
	"""
	def __init__(self, sequence):
		"""Initialize the pipeline with the sequence"""
		self.sequence = sequence

	def generateDeNovoPrediction():
		"""Execute the pipeline and find the minimum conformation"""
		resMapper = ResidueMapper()
		fragLibGen = FragmentLibraryGenerator(self.sequence, resMapper)
		fragLib = fragLibGen.generateLibrary()
		seef = Seef()
		flatChainInit = FlatChainConformationInitializer(self.sequence)
		conformation = flatChainInit.generateConformation()
		conformationSampler = ConformationSampler(conformation, seef, fragLib)
		
		while conformationSampler.hasNext():
			conformation = conformationSampler.next()
			
		min = conformationSampler.minimum()