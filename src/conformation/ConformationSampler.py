from abc import ABCMeta, abstractmethod
import random


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
    def __init__(self, initialConformation, seefModel, fragLib):
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

    """

    def __init__(self, initial_conformation, seef_model, frag_lib):
        """Inits this conformation sampler so iterations can begin"""
        super(ConformationSampler, self).__init__(initial_conformation, seef_model, frag_lib)
        self.conformation = initial_conformation
        self.seef = seef_model
        self.fragLib = frag_lib
        self.minimum_conformation = None
        self.k_max = 5000
        self.k = 0
        self.e_max = -1
        self.e = self.seef.computeEnergy(self.conformation)
        self.temp = 1000
        self.maxTemp = 1000
        self.minTemp = 10
        self.e_best = self.e

    def next_conformation(self):
        """Generates the next conformation using the metropolis algorithm"""
        dummy = self.conformation
        startPos = random.randint(0, len(dummy))

        prob9 = (self.temp - self.minTemp) / (self.maxTemp - self.minTemp)
        count = 9 if prob9 > random.random() else 3
        fragment = self.fragLib.get_kmer_fragments(startPos) if count == 9 else self.fragLib.get_kmer_fragments(startPos)
        for i in range(startPos, startPos + count):
            dummy[i] = fragment[i - startPos]

        energy = self.seef.compute_energy(dummy)

        if P(self.e, energy, self.temp) < random.random():
            self.conformation = dummy
            self.e = energy

        if energy < self.e:
            self.minimum_conformation = dummy
            self.e_best = energy

        self.k += 1
        self.temp -= (self.maxTemp - self.minTemp) / self.k
        return self.conformation

    def has_next(self):
        """Checks if the conditions have been fulfilled for this sampling."""
        return self.k < self.k_max and self.e > self.e_max

    def minimum(self):
        """The current minimum-energy conformation."""
        return self.minimum_conformation
