import sys
from Bio import SeqIO
import subprocess

cdhit_file = sys.argv[1]
cons_file = sys.argv[2]
output = sys.argv[3]

dico = {}
list_seq = []
with open(cdhit_file, "r") as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith(">"):
            try:                
                dico[key] = list_seq
                list_seq = []
            except:            
                pass
        else:            
            if line.endswith("*\n"):
                key = line.split(">")[1].split("...")[0]
            else:
                seq_name = line.split(">")[1].split("...")[0]
                list_seq.append(seq_name)        

seqs = {}
for seq_record in SeqIO.parse(cons_file, "fasta"):
    name = seq_record.id
    seq = seq_record.seq
    seqs[name] = seq


nb_repr_in_cons = 0
nb_singleton = 0
nb_cluster = 0

for seq_id in seqs.keys():
    if seq_id in dico.keys():
        nb_repr_in_cons += 1

        if len(dico[seq_id]) == 0:
            nb_singleton += 1
            with open(output + "/singleton.fa", "a") as f:
                f.write(">" + str(seq_id) + "\n" + str(seqs[seq_id]) + "\n")
        
        else:
            # Work in progress
            nb_cluster += 1            
            subprocess.run(["Refiner"])
            break
                      
print("Number of sequences in the fasta file:\t" + str(len(seqs)))           
print("Number of cluster:\t" + str(len(dico)))
print("Number of representative sequences founded in the fasta file:\t" + str(nb_repr_in_cons))
print("Number of singleton (clusters with only one sequence):\t" + str(nb_singleton))
print("Number of cluster with two or more sequences:\t" + str(nb_cluster))