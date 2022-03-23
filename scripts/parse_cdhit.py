# -*-coding:Utf-8 -*

"""
Author: Tommaso Barberis
Date: 16/03/2022
Description: parse cd-hit-est results, it recover sequences from a cluster and then it calls a consensus using Refiner
"""

import argparse, os
from Bio import SeqIO
import subprocess
import shutil, glob

parser = argparse.ArgumentParser()

parser.add_argument("-i", metavar='', help=('''
cd-hit-est output file.
'''))
parser.add_argument("-fasta", metavar='', help=('''
FASTA file used for the clustering.
'''))
parser.add_argument("-o", default = '.', metavar='', help=('''
output directory for the result file.
'''
))

args = parser.parse_args()

cdhit_file = args.i
copy_file = args.fasta
out_dir = args.o
if out_dir[-1] == "/":
    out_dir = out_dir[:-1]

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def parse_cdhit(cdhit):
    dico = {} # key: representative seq, value: list of the sequences in the cluster
    list_seq = [] # value for the dict

    with open(cdhit, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"): # begin of a new cluster
                try:                
                    dico[key] = list_seq # add the representative sequence as key 
                    list_seq = [] # reset list of sequences
                except:            
                    pass
            else: # within the cluster
                if line.endswith("*\n"): # representative sequence of the cluster
                    key = line.split(">")[1].split("...")[0]
                else: # other sequences
                    seq_name = line.split(">")[1].split("...")[0]
                    list_seq.append(seq_name)        
    return dico

def parse_fasta(fasta):
    seqs = {} # key: name of the copy, value: sequence
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name = seq_record.id
        seq = seq_record.seq
        seqs[name] = str(seq).upper() # it allows to avoid errors while running Refiner
    return seqs

def call_consensus(cdhit_dict, seqs_dict):
    # counters for statistics
    nb_repr_within_copies = 0
    nb_singleton = 0
    nb_cluster = 0

    # output files
    singleton = out_dir + "/singleton.fa"
    consensi = out_dir + "/consensi.fa"
    stk = out_dir + "/consensi.stk"

    c = 0 # counter for naming consensus sequences

    # main
    # if not os.path.isfile(singleton): # if there isn't previously results TODO: uncomment this line and tab the following code block
       
    for seq_id in seqs_dict.keys():
        if seq_id in cdhit_dict.keys():
            nb_repr_within_copies += 1

            if len(cdhit_dict[seq_id]) == 0: # the cluster is singleton
                if not os.path.isfile(singleton): #TODO: delete this line and untab the following block
                    nb_singleton += 1                
                    with open(singleton, "a") as f:
                        f.write(">" + str(seq_id) + "\n" + str(seqs_dict[seq_id]) + "\n")
    
            else: # the cluster is not singleton

                # Work in progress
                nb_cluster += 1
                cluster_file = out_dir + "/" + seq_id + ".clst.fasta"
                seqs_id = cdhit_dict[seq_id]
                seqs_id.append(seq_id) # add representative sequence
                
                # write sequences of the cluster in a fasta file
                with open(cluster_file, "w") as f:
                    for seq in seqs_id:
                        f.write(">" + seq + "\n" + seqs_dict[seq] + "\n")
                
                # call consensus using Refiner
                subprocess.run(["Refiner", cluster_file])
                
                with open(consensi, "a") as f:
                    with open(cluster_file+".refiner_cons", "r") as ref_file:
                        data = ref_file.read()
                        f.write(data)
                
                with open(stk, "a") as f:
                    with open(cluster_file+".refiner.stk", "r") as ref_file:
                        data = ref_file.read()
                        f.write(data)

                # clean the folder
                os.remove(cluster_file + ".refiner_cons")
                os.remove(cluster_file + ".refiner.stk")
                os.remove(cluster_file)
                
                dir_list = glob.iglob(os.path.join(out_dir, "RM_*"))
                for path in dir_list:
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                
                if c == 2: #TODO: remove this block
                    break
                c += 1
    
    with open("stats.txt", "w") as stats:
        stats.write(f"""
        Total number of copies:\t{str(len(seqs_dict))}
        Number of cluster:\t{str(len(cdhit_dict))}
        Number of representative sequences founded in the fasta file with copies:\t{str(nb_repr_within_copies)}
        Number of singleton (clusters with only one sequence):\t{str(nb_singleton)}
        Number of cluster with two or more sequences:\t{str(nb_cluster)}
        """)

if __name__ == "__main__":
    clusters = parse_cdhit(cdhit_file)
    seqs = parse_fasta(copy_file)
    call_consensus(clusters, seqs)
    pass