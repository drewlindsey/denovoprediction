from abc import ABCMeta, abstractmethod
from Residue import *
from ..mapper import map_conformation_to_pdb

def copy_conformation(conformation):
    """Performs a deep copy of the conformation"""
    new_conf = Conformation(conformation.name, [], conformation.experimental)
    for aa in conformation.sequence:
        new_conf.sequence.append(aa)
    
    for residue in conformation.conformation:
        new_conf.conformation.append(residue)
    
    return new_conf

# TODO move conformation initialization to this class ?
class Conformation(object):
    __metaclass__ = ABCMeta
    """A representation of an amino acid tertiary structure.

    Attributes:
        sequence: the amino acid sequence characters
        conformation: A list of Residue values
        pdb_file: The file name
        name: The name of the sequence
    """
    def __init__(self, name, sequence, experimental):
        """Inits Conformation with the conformationInitializer"""
        self.name = name
        self.sequence = sequence
        self.conformation = []
        self.experimental = experimental
        self.pdb_file = None

    def get(self, position):
        """Returns a Residue at the given position"""
        return self.conformation[position]

    def get_length(self):
        """Returns the length of the sequence/conformation/residue list"""
        return len(self.conformation)

    def get_residues(self):
        return self.conformation

    def set(self, position, residue):
        """Sets the angles for this sequence, not changing the amino acid type"""
        residue.set_type(self.sequence[position])
        self.conformation[position] = residue

    def get_name(self):
        """Gets the name of the sequence"""
        return self.name

    @abstractmethod
    def initialize(self):
        """Initializes the backbone for this conformation"""
        pass

    def get_pdb_file(self):
        return self.pdb_file

    def set_pdb_file(self, file_path):
        self.pdb_file = file_path


class LinearBackboneConformation(Conformation):

    def initialize(self):
        """Creates the initial backbone configuration for a linear chain"""
        for i in range(len(self.sequence)):
            self.conformation.append(Residue(self.sequence[i], {"phi": -150, "psi": 150, "omega": 180}))
