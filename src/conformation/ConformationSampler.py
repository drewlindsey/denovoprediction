from __future__ import division

from abc import ABCMeta, abstractmethod
import random
import math
from ..mapper import map_conformation_to_pdb
import copy

#celery worker -A run_webserver.celery --loglevel=info
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
        fragLib: The FragmentLibrary object to sample from
        score_models: The scoring functions, represented by a dictionary {"name1": score_model1, "name2": score_model2}
    """

    @abstractmethod
    def __init__(self, initialConformation, experimentalConformation, seefModel, fragLib, output_loc, score_models):
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

    def __init__(self, initial_conformation, experimental_conformation, seef_model, frag_lib, pdb_output_loc, score_models):
        """Inits this conformation sampler so iterations can begin"""
        super(ConformationSampler, self).__init__(initial_conformation, experimental_conformation, seef_model, frag_lib, pdb_output_loc, score_models)
        self.conformation = initial_conformation
        self.experimental = experimental_conformation
        self.seef = seef_model
        self.scores = score_models
        self.fragLib = frag_lib
        self.minimum_conformation = initial_conformation
        self.k_max = 10
        self.k = 0
        self.e_max = 0
        self.output_loc = pdb_output_loc
        pdb_file = map_conformation_to_pdb(self.conformation, self.output_loc, True)
        self.e = self.seef.compute_score(pdb_file, self.experimental) #self.score.compute_score(pdb_file, self.experimental) #
        self.temp = 2500
        self.maxTemp = 2500
        self.minTemp = 10
        self.e_best = self.e
        self.score_best = 5000
        self.score_best_for_tm = 5000
        self.tm_best_for_score = 0
        self.curr_score = 10000
        print score_models

    def get_k_max(self):
        return self.k_max

    def next_conformation(self):
        """Generates the next conformation using the metropolis algorithm"""
        print self.k
        
        dummy = copy.deepcopy(self.conformation)

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
        #dummy.set_pdb_file(pdb)
        energy = self.seef.compute_score(pdb, self.experimental)
        #energy = self.score.compute_score(pdb, self.experimental)
        self.curr_score = self.scores["rmsd"].compute_score(pdb, self.experimental)
        
        

        # print "[" + str(self.k) + "]" + " ENERGY: " + str(energy)

        probability_acceptance = math.exp(-(self.e - energy) / (self.temp))
        # print "[" + str(self.k) + "]" + " PROB: " + str(probability_acceptance)
        if probability_acceptance > 1 or probability_acceptance > random.random():
            self.conformation = copy.deepcopy(dummy)
            self.conformation.set_pdb_file(pdb)
            self.e = energy

        if energy > self.e_best:
            self.minimum_conformation.set_pdb_file(pdb)
            self.minimum_conformation = copy.deepcopy(dummy)
            self.e_best = energy
            self.score_best_for_tm = self.curr_score
            # print "[" + str(self.k) + "]" + " MINIMUM CHANGE"
            # print "[" + str(self.k) + "]" + " CONFORMATION CHANGE"

        if self.curr_score < self.score_best:
            self.score_best = self.curr_score
            self.tm_best_for_score = energy

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
        output = {}
        for key, value in self.scores:
            output[key] = value.compute_score(map_conformation_to_pdb(self.minimum_conformation, self.output_loc, True), self.experimental)
        return output
        
    def get_current_energy(self):
        """Returns Current Energy"""
        return self.e
        
    def get_current_score(self):
        """Returns Current Score"""
        return self.curr_score
        
    def get_best_energy(self):
    	"""Returns best energy"""
    	return self.e_best
        
    def get_best_score(self):
        """Returns Best Score found"""
        return self.score_best
        
    def get_best_score_for_tm(self):
        """Resturns Score For Best Tm-Score"""
        return self.score_best_for_tm
        
    def get_best_tm_for_score(self):
        """Returns Tm-score for best score"""
        return self.tm_best_for_score
