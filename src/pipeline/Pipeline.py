from abc import ABCMeta, abstractmethod
from ..fragments.FragmentLibrary import RobettaFragmentLibrary
from ..seef.Seef import *
from ..score.Score import *
from ..conformation.ConformationSampler import *
from ..conformation.Conformation import *


class BasePipeline(object):
    __metaclass__ = ABCMeta
    """An abstract base class for a de novo pipeline. An implementation should have all
    aspects of the pipeline to generate the predicted structure.

    TODO: pass in the required objects to build this pipeline??

    Attributes:
        sequence: The amino acid sequence
    """

    @abstractmethod
    def __init__(self, sequence):
        """Inits the pipeline with the given sequence"""
        self.conformation = None
        self.experimental = None
        self.complete = False
        self.result = 0
        pass

    @abstractmethod
    def generate_structure_prediction(self, output_loc):
        """Execute the pipeline and obtain the predicted structure conformation"""
        pass

    def get_current_conformation(self):
        return self.conformation

    def is_complete(self):
        return self.complete


class LinearPipeline(BasePipeline):
    """A linear backbone initialized default pipeline.

    Attributes:
        sequence: The sequence to determine the tertiary structure
        robetta_dict: a dictionary containing the k (for k-mer) and the associated
        robetta text file for that k-value
    """

    def __init__(self, name, sequence, robetta_dict, experimental_pdb):
        """Initialize the pipeline with the sequence"""
        super(LinearPipeline, self).__init__(sequence)
        self.name = name
        self.sequence = sequence
        self.robetta_dict = robetta_dict
        self.experimental = experimental_pdb

    def generate_structure_prediction(self, output_loc):
        """Execute the pipeline and find the minimum conformation"""
        frag_lib = RobettaFragmentLibrary(self.sequence)
        frag_lib.generate(self.robetta_dict)
        seef = DFirePotential()
        score = TMScore()
        self.conformation = LinearBackboneConformation(self.name, self.sequence)
        self.conformation.initialize()
        sampler = ConformationSampler(self.conformation, self.experimental, seef, score, frag_lib, output_loc)

        while sampler.has_next():
            self.conformation = sampler.next_conformation()
            # TODO submit to 3dmol.js

        self.conformation = sampler.minimum()
        self.complete = True
        self.result = sampler.score_conformation()
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n" \
              "%%%           PIPELINE COMPLETE           %%%\n" \
              "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
