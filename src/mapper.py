def map_robetta_structure_to_fragments(k, input_file):
    """Takes the input file and outputs the k-mer fragment list"""
    fragments = []
    start = False
    count = 0
    frag_count = 0
    with open(input_file) as robFile:
        for line in robFile:
            if count == 2 and start:
                start = True
                count = 0
            if count == 4 and not start:
                count = 0
                if frag_count == 200:
                    frag_count = 0
            if count < 4 and not start:
                frag_count += 1
                data = line[18:44]
                print data
            count += 1

    return fragments


def map_conformation_to_pdb(conformation):
    pass
