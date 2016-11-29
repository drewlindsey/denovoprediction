import re
import subprocess as sp
import os
import tempfile
import uuid
from fragments.Fragment import Fragment
from conformation.Residue import Residue


def map_robetta_structure_to_fragments(k, input_file):
    """Takes the input file and outputs the k-mer fragment list

    The k-mer fragment list appears like this (n = length of sequence - k, position = 200)
    Each frag_ij for index i and position j is a k-mer fragment

                index
                0           1           2         ...           n

    position  0 frag_00     frag_01     frag_02                 frag_0n
              1 frag_10     frag_11     frag_12                 frag_1n
              2 frag_20     frag_21     frag_22                 frag_2n
              3 frag_30     frag_31     frag_32                 frag_3n
              ...
              m frag_m0     frag_m1     frag_m2                 frag_mn

    The algorithm works by iterating over the 200 k-mer fragments for each position adding each k-mer
    fragment to the list at the position
    """
    fragments = []
    count = 0
    frag_count = 0
    curr_position = 0
    k_list = []
    fragments.append([])
    with open(input_file) as robFile:
        for line in robFile:
            if frag_count == 200:
                fragments.append([])
                frag_count = 0
                curr_position += 1

            if count >= k:
                fragment = Fragment(k_list, k)
                fragments[curr_position].append(fragment)
                k_list = []
                count = 0
                frag_count += 1
            angle_line = line[18:44]
            result = re.findall('(-?[0-9]+.[0-9]+)', angle_line)

            # we found a match to the regex
            if len(result) > 0:
                if len(result) != 3:
                    continue
                residue = Residue(None, {'phi': result[0], 'psi': result[1], 'omega': result[2]})
                k_list.append(residue)
                count += 1

    return fragments


def map_conformation_to_pdb(conformation, loc):
    """Takes a conformation (Conformation.py) object and creates a PDB file (using crankite)
    and return the path to this file.

    Arguments:
        conformation: the Conformation object
        loc: the directory to store the pdb file
    """

    tmp_file_name = loc + conformation.name + "-" + str(uuid.uuid1()) + "_tmp.txt"
    with open(tmp_file_name, 'w+') as tmp_file:
        for residue in conformation.get_residues():
            angles = residue.get_angles()
            tmp_file.write("{0} {1} {2} {3}\n".
                           format(residue.get_type(), angles["phi"], angles["psi"], angles["omega"]))

    print "IS FILE ??? " + str(os.path.isfile(tmp_file_name))

    file_name = os.path.join(loc, conformation.get_name() + "-" + str(uuid.uuid1()) + ".pdb")

    with open(file_name) as some_file:
        lipa_call = sp.Popen(['lipa', tmp_file_name, '>', file_name], stdout=sp.PIPE, stderr=sp.PIPE)
        lipa_out, err = lipa_call.communicate()

        print lipa_out

    return file_name
