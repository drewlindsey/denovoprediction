from abc import ABCMeta, abstractmethod
import random
import math
from ..mapper import map_conformation_to_pdb


class BaseConformationSampler(object):
    __metaclass__ = ABCMeta
    """An abstract base class for the ConformationSampler.

    Controls sampling the conformation space and testing the energy. Acts
    as though it is an iterator and can be continually iterated upon until
    the hasNext conditions have been completed.

    Attributes:
        initialConformation: The initial backbone conformation
        seefModel: The SEEF
        fragLib: The FragmentLibrary object to sample from
    """

    @abstractmethod
    def __init__(self, initialConformation, seefModel, fragLib, output_loc):
        """Inits the ConformationSampler"""
        pass

    @abstractmethod
    def next_conformation(self):
        """Generate the next conformation"""
        pass

    @abstractmethod
    def has_next(self):
        """Checks if the sampling is complete"""
        pass

    @abstractmethod
    def minimum(self):
        """The current minima conformation"""
        pass


class ConformationSampler(BaseConformationSampler):
    """An implementation of AbstractConformationSampler. Uses the SEEF and fragment library to generate a new sample and test
    its energy. End goal is to determine the minimum energy conformation. Works similar to an iterator.

    Attributes:
        initialConformation: The initial backbone conformation
        seefModel: The SEEF
        fragLib: The k-mer neighbor fragment library
        pdb_output_loc: location to store pdb files
    """

    def __init__(self, initial_conformation, seef_model, frag_lib, pdb_output_loc):
        """Inits this conformation sampler so iterations can begin"""
        super(ConformationSampler, self).__init__(initial_conformation, seef_model, frag_lib, pdb_output_loc)
        self.conformation = initial_conformation
        self.seef = seef_model
        self.fragLib = frag_lib
        self.minimum_conformation = initial_conformation
        self.k_max = 10
        self.k = 1
        self.e_max = -20000
        self.output_loc = pdb_output_loc
        self.e = self.seef.compute_energy(map_conformation_to_pdb(self.conformation, self.output_loc))
        self.temp = 1000
        self.maxTemp = 1000
        self.minTemp = 10
        self.e_best = self.e

    def next_conformation(self):
        """Generates the next conformation using the metropolis algorithm"""
        dummy = self.conformation
        print "dummy length : " + str(dummy.get_length())

        prob9 = (self.temp - self.minTemp) / (self.maxTemp - self.minTemp)
        count = 9 if prob9 > random.random() else 3

        startPos = random.randint(0, dummy.get_length()-(count+1))
        rand_neighbor = random.randint(0, 199)

        print "start: " + str(startPos)
        print "count: " + str(count)
        print "length of kmer fragment list " + str(len(self.fragLib.get_kmer_fragments(count, startPos)))

        # gets a residue in the (startPos,rand_neighbor position)
        fragment = self.fragLib.get_kmer_fragment(count, startPos, rand_neighbor)

        for i in range(startPos, startPos + count):
            print "i: " + str(i)
            print "i-start " + str(i-startPos)
            # assign the residue
            dummy.set(i, fragment.get_residue(i - startPos))

        energy = self.seef.compute_energy(map_conformation_to_pdb(dummy, self.output_loc))
        print "[" + str(self.k) + "]" + " ENERGY: " + str(energy)

        probability_acceptance = math.exp(-(self.e - energy)) / (self.k * self.temp)
        if probability_acceptance > random.random():
            print "[" + str(self.k) + "]" + " CONFORMATION CHANGE"
            self.conformation = dummy
            self.e = energy

        if energy < self.e:
            print "[" + str(self.k) + "]" + " MINIMUM CHANGE"
            self.minimum_conformation = dummy
            self.e_best = energy

        self.k += 1
        print "[" + str(self.k) + "]" + " TEMP: " + str(self.temp)
        self.temp -= (self.maxTemp - self.minTemp) / self.k_max

        return self.conformation

    def has_next(self):
        """Checks if the conditions have been fulfilled for this sampling."""
        print "[" + str(self.k) + "]" + " HAS NEXT: " + str(self.k < self.k_max and self.e > self.e_max)
        return self.k < self.k_max and self.e > self.e_max

    def minimum(self):
        """The current minimum-energy conformation."""
        return self.minimum_conformation
