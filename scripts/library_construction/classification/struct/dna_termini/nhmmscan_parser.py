"""
Author: Tommaso BARBERIS
Date: 24/06/2022
Description: script to parse results from nhmmscan (obtained using profiles of DNA 
    termini provided by Dfam: https://www.dfam.org/releases/dna_termini_1.1/dna_termini_1.1.hmm.gz)

Usage: python3 nhmmscan_parser.py -i file.tsv 
Used protocol:
hmmpress dna_termini_1.1.hmm
nhmmscan -E 1e-2 --tblout hmmscan.tsv dna_termini_1.1.hmm consensi_rnd8.renamed.polished.fa

For python dependancies you can use the conda env TE_Aalb (https://github.com/TommasoBarberis/TE_Aalb)
"""

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", metavar='', help=('''TSV file obtained from nhmmscan'''))
    

    args = parser.parse_args()

    # assign parameters
    tsv_file = args.i

    hits_dict = {}
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#"):
                pass
            else:
                fields = line.split()
                target = fields[0] # it correspond to the possible classification
                seq = fields[2]

                if seq not in hits_dict.keys():
                    hits_dict[seq] = [target]
                else:
                    tmp_lst = hits_dict[seq]
                    tmp_lst.append(target)
                    hits_dict[seq] = tmp_lst


    with open("dna_classified", "w") as f:
        for seq in hits_dict.keys():
            hits = hits_dict[seq]
            hits = [x.replace("_begins", "").replace("_ends", "").replace("_termini", "") for x in hits]

            if len(set(hits)) > 1:
                # handle sequences with differents hits
                hits = [x.split("_") for x in hits]
                intersection = set.intersection(*map(set,hits)) # remove the divergent part of the classification
                    
                order_dict = {}
                for val in intersection: # it allows to call the differents terms in the right order
                    ind = hits[0].index(val)
                    order_dict[ind] = val
    
                c = 0
                classif = []
                while True:
                    if c in order_dict.keys():
                        classif.append(order_dict[c])
                    else:
                        break
                    c += 1
                if len(classif) > 0:
                    classif = "-".join(classif)
                else:
                    classif = "Unknown"
                
            else:
                # no ambigous hits
                classif = hits[0].replace("_", "-")
            
            f.write(classif + "\t" + seq + "\n")