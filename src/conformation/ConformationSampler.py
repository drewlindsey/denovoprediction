from __future__ import division

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
        experimentalConformation: The experimental native structure that we are trying to predict
        seefModel: The SEEF
        scoreModel: The scoring function
        fragLib: The FragmentLibrary object to sample from
    """

    @abstractmethod
    def __init__(self, initialConformation, experimentalConformation, seefModel, scoreModel, fragLib, output_loc):
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
        
    @abstractmethod
    def score_conformation(self):
        """Score the conformation"""
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

    def __init__(self, initial_conformation, experimental_conformation, seef_model, score_model, frag_lib, pdb_output_loc):
        """Inits this conformation sampler so iterations can begin"""
        super(ConformationSampler, self).__init__(initial_conformation, experimental_conformation, seef_model, score_model, frag_lib, pdb_output_loc)
        self.conformation = initial_conformation
        self.experimental = experimental_conformation
        self.seef = seef_model
        self.score = score_model
        self.fragLib = frag_lib
        self.minimum_conformation = initial_conformation
        self.k_max = 5000
        self.k = 0
        self.e_max = -20000
        self.output_loc = pdb_output_loc
        pdb_file = map_conformation_to_pdb(self.conformation, self.output_loc, True)
        self.e = self.seef.compute_energy(pdb_file)
        self.temp = 5000
        self.maxTemp = 5000
        self.minTemp = 10
        self.e_best = self.e

    def get_k_max(self):
        return self.k_max

    def next_conformation(self):
        """Generates the next conformation using the metropolis algorithm"""
        print self.k
        
        dummy = self.conformation

        prob9 = (self.temp - self.minTemp) / (self.maxTemp - self.minTemp)
        count = 9 if prob9 > random.random() else 3

        self.temp -= (self.maxTemp - self.minTemp) / self.k_max
        # print "[" + str(self.k) + "]" + " TEMP: " + str(self.temp)

        startPos = random.randint(0, dummy.get_length() - (count + 1))
        rand_neighbor = random.randint(0, 199)

        # gets a residue in the (startPos,rand_neighbor position)
        fragment = self.fragLib.get_kmer_fragment(count, startPos, rand_neighbor)

        for i in range(startPos, startPos + count):
            # assign the residue
            dummy.set(i, fragment.get_residue(i - startPos))

        pdb = map_conformation_to_pdb(dummy, self.output_loc, True)
        dummy.set_pdb_file(pdb)
        energy = self.seef.compute_energy(pdb)

        # print "[" + str(self.k) + "]" + " ENERGY: " + str(energy)

        probability_acceptance = math.exp(-(self.e - energy) / (self.temp))
        # print "[" + str(self.k) + "]" + " PROB: " + str(probability_acceptance)
        if probability_acceptance > 1 or probability_acceptance > random.random():
            self.conformation = dummy
            self.e = energy

        if energy < self.e_best:
            self.minimum_conformation = dummy
            self.e_best = energy
            # print "[" + str(self.k) + "]" + " MINIMUM CHANGE"
            # print "[" + str(self.k) + "]" + " CONFORMATION CHANGE"

        # print "[" + str(self.k) + "]" + " BEST ENERGY SO FAR: " + str(self.e_best)
        self.k += 1

        return self.conformation

    def has_next(self):
        """Checks if the conditions have been fulfilled for this sampling."""
        # print "[" + str(self.k) + "]" + " HAS NEXT: " + str(self.k < self.k_max and self.e > self.e_max)
        return self.k < self.k_max and self.e > self.e_max

    def minimum(self):
        """The current minimum-energy conformation."""
        return self.minimum_conformation
        
    def score_conformation(self):
        """Score the conformation"""
        sc = self.score.compute_score(map_conformation_to_pdb(self.minimum_conformation, self.output_loc, True), self.experimental)
        #print "SCORE: " + str(sc)
        return sc        
