from abc import ABCMeta, abstractmethod


# TODO move conformation initialization to this class ?
class Conformation(object):
    __metaclass__ = ABCMeta
    """A representation of an amino acid tertiary structure.

    Attributes:
        sequence: the amino acid sequence characters
        conformationInitializer: An object of type AbstractConformationInitializer to create the initial "backbone" (list of Residues)
        conformation: A list of Residue values
    """
    def __init__(self, sequence):
        """Inits Conformation with the conformationInitializer"""
        self.sequence = sequence
        self.conformation = self.conformation_initializer.generateConformation()

    def reset(self, initializer=None):
        """Resets back to the initial conformation using the original ConformationInitializer"""
        if initializer is not None:
            self.conformation_initializer = initializer
        self.conformation = self.conformation_initializer.generateConformation()

    def get(self, position):
        """Returns a Residue at the given position"""
        return self.conformation[position]

    @abstractmethod
    def initialize(self):
        pass


class LinearBackboneConformation(Conformation):
    def initialize(self):
        """Creates the initial backbone configuration for a linear chain"""
        pass
