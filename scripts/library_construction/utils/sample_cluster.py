"""
Author: Tommaso Barberis
Date: 08/04/2022
Description: from the output of cdhit, it samples some clusters
"""

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from random import sample
import argparse
from parse_cdhit import parse_cdhit, parse_fasta


parser = argparse.ArgumentParser()

parser.add_argument("--cdhit", metavar='', help=('''
cd-hit-est output file.
'''))
parser.add_argument("--fasta", metavar='', help=('''
FASTA file with copies used for the clustering with cd-hit-est.
'''))
parser.add_argument("-n", metavar='', default='10', help=('''
number of cluster to sample.
'''))
parser.add_argument("--min_size", metavar='', default='5', help=('''
minimum size of the cluster.
'''))
parser.add_argument("--max_size", metavar='', default='0', help=('''
maximum size of the cluster, when equal 0, no limit is set.
'''))

args = parser.parse_args()

cdhit_file = args.cdhit
fasta_file = args.fasta
nb_cluster = int(args.n)
min_size = int(args.min_size)
max_size = int(args.max_size)

def filter(cdhit_file, min_size, max_size):
    """
    From a cdhit output file return a dictionary with clusters that have a minimum and a maximum size.
    """  
    cdhit_dict = parse_cdhit(cdhit_file)
    new_dict = {}

    for clst in cdhit_dict.keys():
        clst_size = len(cdhit_dict[clst])
        if clst_size > min_size:
            if max_size == 0 or clst_size < max_size:
                new_dict[clst] = cdhit_dict[clst]
            else:
                pass
        else:
            pass
    return new_dict

def create_samples(clst_sampled, cdhit_dict, seq_dict):
    """
    For the N clusters it creates samples of size: 
    100, 150, 200, 250,300, 350, 400, 450, 500
    and then it creates also samples with 500 sequences with a given lenght threshold (by the maximum length):
    0.2, 0.4, 0.6, 0.8, 1
    """

    c = 0 # counter, to name folders
    sizes_lst = [100, 150, 200, 250,300, 350, 400, 450, 500] # differentes subsamples in the cluster
    
    for clst in clst_sampled:
        
        seq_id_lst = cdhit_dict[clst]
        seq_id_lst.append(clst)

        # create the folder 
        os.mkdir("./cluster_" + str(c))

        cluster_file_prefix = "./cluster_" + str(c) + "/"

        # for all sample sizes
        for s in sizes_lst:
            # the file name
            cluster_file = cluster_file_prefix + "cluster-" + str(s) + ".clst.fasta"
            counter_seq = 0

            for seq_id in seq_id_lst:
            
                with open(cluster_file, "a+") as f:        
                    f.write(">" + "sequence_" + str(counter_seq) + "\n" + seq_dict[seq_id] + "\n")
                    counter_seq += 1

                    # if the number of sample size is reach, close the file
                    if counter_seq == s:
                        break

            if s == 500:
                max_len = 0
                len_lst = [2, 4, 6, 8, 10]
                
                for seq_id in seq_id_lst:
                    seq_len = len(seq_dict[seq_id])
                    
                    if seq_len > max_len:
                        max_len = seq_len

                for l in len_lst:
                    cluster_file = cluster_file_prefix + "cluster-" + str(l) + ".clst.fasta"

                    for seq_id in seq_id_lst:
                        with open(cluster_file, "a+") as f:     
                            seq_len = len(seq_dict[seq_id])
                            threshold = l * max_len / 10

                            if seq_len >= threshold:
                                f.write(">" + "sequence_" + str(counter_seq) + "\n" + seq_dict[seq_id] + "\n")
                                counter_seq += 1

                                # if the number of sample size is reach, close the file
                                if counter_seq == s:
                                    break

        c += 1



if __name__ == "__main__":
    filter_dict = filter(cdhit_file, min_size, max_size)
    clst_sampled = sample(list(filter_dict.keys()), nb_cluster)
    seq_dict = parse_fasta(fasta_file)
    create_samples(clst_sampled, filter_dict, seq_dict)
         