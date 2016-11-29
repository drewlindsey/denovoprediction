from abc import ABCMeta, abstractmethod
import subprocess as sp
import re

class BaseScore(object):
	__metaclass__ = ABCMeta
	"""An abstract class template for the RMSD and TM_Score Functions
	
	Attributes:
		options: the options for the score functions (if any)
	"""
	
	@abstractmethod
	def __init__(self, *args):
		"""Inits the Score with the given (optional) options"""
		pass
		
	@abstractmethod
	def compute_score(self, conformation, experimental):
		"""Computes the score of the given conformation vs. an 
		experimental native structure"""
		pass
		
class TMScore(BaseScore):
	"""A Scoring function requiring linux and the TMscore 
	(http://zhanglab.ccmb.med.umich.edu/TM-score/) executable
	to be downloaded and compiled using the fortran compiler
	with the command '>gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f'
	on the PATH or in the directory calling this function.
	"""
	
	def __init__(self):
		super(TMScore, self).__init__()
		
	def compute_score(self, conformation, experimental):
		"""Generates the TMScore for the given conformation using TMscore executable file
		
		conformation:
			A PDB file to compare with an experimental native structure
			
		experimental:
			A PDB file with to compare the given conformation against
		"""
		
		tmScoreCall = sp.Popen(['./TMscore', conformation, experimental], stdout = sp.PIPE, stderr = sp.PIPE)
		tmScoreString, err = tmScoreCall.communicate()
		tmScoreValueString = re.search("TM-score\s+=\s([0-9]+\.[0-9]+)", tmScoreString).group(0)
		tmScoreValue = re.findall("\d+\.\d+", tmScoreValueString)
		return float(tmScoreValue[0])
		
class RMSDScore(BaseScore):
	"""A Scoring function requiring linux and the TMScore 
	(http://zhanglab.ccmb.med.umich.edu/TM-score/) executable
	to be downloaded and compiled using the fortran comipler
	with the command '>gfortran -static -O3 -ffast-math -lm -o TMscore TMscore.f'
	on the PATH or in the directory calling this function.
	"""
	
	def __init__(self):
		super(RMSDScore, self).__init__()
		
	def computer_score(self, conformation, experimental):
		"""Generates the RMSDScore for the given conformation using TMscore executable file
		
		conformation:
			A PDB file to compare with an experimental native structure
			
		experimental:
			A PDB file to compare the given conformation against
		"""
		
		rmsdScoreCall = sp.Popen(['./TMscore', conformation, experimental], stdout = sp.PIPE, stderr = sp.PIPE)
		rmsdScoreString, err = rmsdScoreCall.communicate()
		#RMSD of  the common residues=    0.000
		rmsdScoreValueString = re.search("RMSD\s+of\s+the\s+common\s+residues=\s+([0-9]+\.[0-9]+)", rmsdScoreString).group(0)
		rmsdScoreValue = re.findall("\d+\.\d+", rmsdScoreValueString)
		return float(rmsdScoreValue[0])
		
		
		
		
		
		
		
		
		
		
		
