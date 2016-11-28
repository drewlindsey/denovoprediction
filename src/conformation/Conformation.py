from abc import ABCMeta, abstractmethod
from Residue import *
from ..mapper import map_conformation_to_pdb


# TODO move conformation initialization to this class ?
class Conformation(object):
    __metaclass__ = ABCMeta
    """A representation of an amino acid tertiary structure.

    Attributes:
        sequence: the amino acid sequence characters
        conformation: A list of Residue values
        pdb_mapper: The PDB mapper object
        pdb_file: The raw text PDB file for this conformation
    """
    def __init__(self, sequence, pdb_mapper):
        """Inits Conformation with the conformationInitializer"""
        self.sequence = sequence
        self.conformation = []
        self.pdb_mapper = pdb_mapper
        self.pdb_file = None

    def get(self, position):
        """Returns a Residue at the given position"""
        return self.conformation[position]

    @abstractmethod
    def initialize(self):
        """Initializes the backbone for this conformation"""
        pass

    @abstractmethod
    def get_pdb_file(self):
        """Updates and gets the pdb file for the current conformation"""


class LinearBackboneConformation(Conformation):
    def get_pdb_file(self):
        self.pdb_file = map_conformation_to_pdb(self.conformation)

    def initialize(self):
        """Creates the initial backbone configuration for a linear chain"""
        for i in range(len(self.sequence)):
            self.conformation[i] = Residue(self.sequence[i], {"phi": 180, "psi": -180, "omega": 180})
