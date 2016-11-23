from abc import ABCMeta, abstractmethod


class BaseSeef(object):
	__metaclass__ = ABCMeta
	"""An abstract class template for the Statistical Effective Energy Function (SEEF)
	
	Attributes:
		options: the options for the SEEF function (if any)
	"""

	@abstractmethod
	def __init__(self, options=None):
		"""Inits the SEEF with the given (optional) options"""
		pass

	@abstractmethod
	def compute_energy(self, conformation):
		"""Computes the energy of the given conformation as described by the implementation details"""
		pass
