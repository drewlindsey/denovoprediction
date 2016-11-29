from abc import ABCMeta, abstractmethod
from ..mapper import *


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
        self.fragments = {}

    def get_kmer_fragments(self, k, index):
        """Gets the k-mer Fragment 1D list for the given sequence index"""
        return self.fragments[k][index]

    def get_kmer_fragment(self, k, index, position):
        """Gets a k-mer Fragment for the given sequence index at the given position in the fragment list"""
        print k
        print index
        print len(self.fragments[k])
        print position
        print len(self.fragments[k][index])

        return self.fragments[k][index][position]

    @abstractmethod
    def generate(self, file_dict):
        """Generate a fragment library from the given file_dict
        which should be {k : file} pairs"""
        pass


class RobettaFragmentLibrary(BaseFragmentLibrary):
    """Generates a fragment library for a robetta fragment"""

    def generate(self, file_dict):
        """Generates fragments using the robetta mapper."""
        for key in file_dict:
            self.fragments[int(key)] = map_robetta_structure_to_fragments(key, file_dict[key])
