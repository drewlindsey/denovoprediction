from abc import ABCMeta, abstractmethod
from ..fragments.FragmentLibrary import BaseFragmentLibrary
from ..seef.Seef import *
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
        pass

    @abstractmethod
    def generate_structure_prediction(self):
        """Execute the pipeline and obtain the predicted structure conformation"""
        pass


class LinearPipeline(BasePipeline):
    """A linear backbone initialized default pipeline.

    Attributes:
        sequence: The sequence to determine the tertiary structure
    """

    def __init__(self, sequence):
        """Initialize the pipeline with the sequence"""
        super(LinearPipeline, self).__init__(sequence)
        self.sequence = sequence

    def generate_structure_prediction(self):
        """Execute the pipeline and find the minimum conformation"""
        frag_lib = BaseFragmentLibrary(self.sequence)
        seef = BaseSeef()
        conformation = LinearBackboneConformation(self.sequence)
        conformation.initialize()
        sampler = ConformationSampler(conformation, seef, frag_lib)

        while sampler.hasNext():
            conformation = sampler.next_conformation()
            # TODO submit to 3dmol.js

        min_conf = sampler.minimum()
