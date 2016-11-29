from abc import ABCMeta, abstractmethod
import subprocess as sp
import re
import os


class BaseSeef(object):
    __metaclass__ = ABCMeta
    """An abstract class template for the Statistical Effective Energy Function (SEEF)

    Attributes:
        options: the options for the SEEF function (if any)
    """

    @abstractmethod
    def __init__(self, *args):
        """Inits the SEEF with the given (optional) options"""
        pass

    @abstractmethod
    def compute_energy(self, conformation):
        """Computes the energy of the given conformation as described by the implementation details"""
        pass


class RwPotential(BaseSeef):
    """A SEEF function requiring linux and the calRw (http://zhanglab.ccmb.med.umich.edu/RW/)
    executable to be installed on the PATH or in the directory calling this function.
    """

    def __init__(self):
        super(RwPotential, self).__init__()

    def compute_energy(self, conformation_pdb):
        """Generates the energy value for the given conformation using RW Potential executable file

        conformation_pdb:
            A PDB file to determine energy from
        """

        rwPotentialCall = sp.Popen(['./calRW', conformation_pdb], stdout=sp.PIPE, stderr=sp.PIPE)
        rwPotentialString, err = rwPotentialCall.communicate()
        rwPotentialValue = re.search("-?([0-9]+\.[0-9]+)", rwPotentialString)
        if rwPotentialValue is None:
            return float("inf")
        return float(rwPotentialValue.group(0))


class DFirePotential(BaseSeef):
    """A SEEF function requiring linux and the dDFIRE (http://sparks-lab.org/yueyang/DFIRE/dDFIRE-service.php)
    executable to be installed on the PATH or in the directory calling this function.
    """

    def __init__(self):
        super(DFirePotential, self).__init__()

    def compute_energy(self, conformation):
        """Generates the energy value for the given conformation using the dDFIRE executable file

        conformation:
            A PDB file to determine engery from
        """
        dDFireCall = sp.Popen(['./dDFIRE', conformation], stdout=sp.PIPE, stderr=sp.PIPE)
        dDFireString, err = dDFireCall.communicate()
        dDFireValue = re.search("(-?[0-9]+\.[0-9]+)", dDFireString).group(0)
        return float(dDFireValue)
