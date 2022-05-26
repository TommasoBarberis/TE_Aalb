"""
Author: Tommaso BARBERIS
Date: 25/05/22
Description: parse RepeatMasker .out file and return tsv file with the best hit for each query sequence
"""

import argparse


def parse_rmout(rm_file):
    
    hits = {}
    with open(rm_file, 'r') as f:
        data = f.readlines()
        for line in data:
            if line.startswith("   SW") or line.startswith("score") or len(line) == 1:
                pass
            else:
                fields = line.split()
                if fields[4] not in hits.keys():
                    hits[fields[4]] = (fields[0], fields[1], fields[9], fields[10])
                else:
                    old_best_hit = hits[fields[4]]
                    if fields[0] > old_best_hit[0]:
                        hits[fields[4]] = (fields[0], fields[1], fields[9], fields[10])
    
    with open(rm_file + ".tsv", "w") as f:
        f.write("query_name\tscore_SW\tdivergence\tmatching_repeat\trepeat_class_family\n")
        for seq in hits.keys():
            f.write(seq + "\t" + hits[seq][0] + "\t" + hits[seq][1] + "\t" + hits[seq][2] + "\t" + hits[seq][3] + "\n")
                


if  __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-r", metavar='', help=('''
    RepeatMasker .out file
    '''))
    
    args = parser.parse_args()

    # assign parameters
    rm_file = args.r
    parse_rmout(rm_file)
    pass