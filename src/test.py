import os
import sys

cwd = os.getcwd()
#sys.path.insert(0, cwd+'\\mapping')

from mapping.FragmentMapper import FragmentMapper

mapper = FragmentMapper(3)
mapper.mapFileToFragments('C:\\Users\\clindsd8\\Downloads\\CASP11_T0856_HERC1\\CASP11_T0856_HERC1\\aat000_03_05.200_v1_3.txt')
