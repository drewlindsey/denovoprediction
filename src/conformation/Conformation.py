from abc import ABCMeta, abstractmethod
from Residue import *


# TODO move conformation initialization to this class ?
class Conformation(object):
    __metaclass__ = ABCMeta
    """A representation of an amino acid tertiary structure.

    Attributes:
        sequence: the amino acid sequence characters
        conformation: A list of Residue values
    """
    def __init__(self, sequence):
        """Inits Conformation with the conformationInitializer"""
        self.sequence = sequence
        self.conformation = []

    def get(self, position):
        """Returns a Residue at the given position"""
        return self.conformation[position]

    @abstractmethod
    def initialize(self):
        """Initializes the backbone for this conformation"""
        pass


class LinearBackboneConformation(Conformation):
    def initialize(self):
        """Creates the initial backbone configuration for a linear chain"""
        for i in range(len(self.sequence)):
            self.conformation[i] = Residue(self.sequence[i], {"phi": 180, "theta": -180, "omega": 180})
