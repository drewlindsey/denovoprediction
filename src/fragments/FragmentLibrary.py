from abc import ABCMeta, abstractmethod


class BaseFragmentLibrary(object):
    __metaclass__ = ABCMeta
    """A representation of a library of fragments

    This library contains both 3-mers and 9-mers which are lists containing at each index
    of the sequence a 3 or 9 residue fragment for that position.

    Attributes:
        sequence: The target sequence
        frag3: The 3-mer fragment lists (2-D array of Fragment objects) (dimensions are typically length x 200)
        frag9: The 9-mer fragment lists (2-D array of Fragment objects) (dimensions are typically length x 200)
    """

    def __init__(self, sequence):
        """Initializes the FragmentLibrary for a sequence of length length with the 2D fragment lists"""
        self.sequence = sequence
        self.fragments3 = []
        self.fragments9 = []

    def get_3mer_fragments(self, index):
        """Gets the 3-mer Fragment 1D list for the given sequence index"""
        return self.fragments3[index]

    def get_9mer_fragments(self, index):
        """Gets the 9-mer Fragment 1D list for the given sequence index"""
        return self.fragments9[index]

    def get_3mer_fragment(self, index, position):
        """Gets a 3-mer Fragment for the given sequence index at the given position in the fragment list"""
        return self.fragments3[index][position]

    def get_9mer_fragment(self, index, position):
        """Gets a 9-mer Fragment for the given sequence index at the given position in the fragment list"""
        return self.fragments9[index][position]

    @abstractmethod
    def generate(self, file_name, residue_mapper):
        pass


class RobettaFragmentLibrary(BaseFragmentLibrary):
    def generate(self, file_name, residue_mapper):
        pass
