class LinearPipeline(AbstractPipeline):
	def __init__(self, sequence):
		self.sequence = sequence

	def generateDeNovoPrediction():
		resMapper = ResidueMapper()
		fragLibGen = FragmentLibraryGenerator(self.sequence, resMapper)
		fragLib = fragLibGen.generateLibrary()
		seef = Seef()
		flatChainInit = FlatChainConformationInitializer(self.sequence)
		conformation = flatChainInit.generateConformation()
		conformationSampler = ConformationSampler(conformation, seef, fragLib)
		
		while(conformationSampler.hasNext())
			conformation = conformationSampler.next()
			
		min = conformationSampler.minimum()