from abc import ABCMeta, abstractmethod


class BaseConformationInitializer(object):
    __metaclass__ = ABCMeta
    """An abstract base class for conformation initialization. Determines
    how the conformation is initialized, e.g. straight-line, 'nearly' straight-line,
    or a 'match' from the PDB

    Attributes:
        sequence: The sequence to initialize into a list of Residue objects
    """

    @abstractmethod
    def __init__(self, sequence):
        """Inits the Conformation Initializer"""
        pass

    def generate_conformation(self):
        """Generates the initialized conformation using this class' methods"""
        pass


class FlatChainConformationInitializer(BaseConformationInitializer):
    """Initializes a conformation with linear, 180 degree, angles

    Attributes:
        sequence: The sequence to initialize
    """

    def __init__(self, sequence):
        """Inits the conformation initializer."""
        super(FlatChainConformationInitializer, self).__init__(sequence)
        self.sequence = sequence

    def generate_conformation(self):
        """Generates a straight-line backbone conformation."""
        raise Exception('Not implemented yet, cuz lazy')
