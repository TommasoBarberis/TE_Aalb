# -*-coding:Utf-8 -*

"""
Author: Tommaso Barberis

Date: 16/03/2022

Description: parse cd-hit-est results. It recovers FASTA sequences used to call cd-hit-est for each cluster. Then it 
separates singleton clusters in 'singleton.fa', low copies number families in 'low_copy_number.fa' and it calls a consensus 
for all other clusters using Refiner (from RepeatModeler2). Consensi sequences can founded in 'consensi.fa' and alignements
in 'consensi.stk'.
The 'stats.txt' reports general statitics on clusters.

Warning: temporary file from Refiner are not conserved.

Dependancies:
    - Refiner 

Tips:
You can install python dependencies using the .yml in the root directory of the project:
conda env create -f TE_Aalb.yml
And then make them available using:
conda activate TE_Aalb
"""

import argparse, os
from Bio import SeqIO
import subprocess, time
import shutil, glob
import multiprocessing as mp
import plotly.express as px

parser = argparse.ArgumentParser()

parser.add_argument("--cdhit", metavar='', help=('''
cd-hit-est output file.
'''))
parser.add_argument("--fasta", metavar='', help=('''
FASTA file used for the clustering.
'''))
parser.add_argument("-o", default = '.', metavar='', help=('''
output directory for the result file.
'''
))
parser.add_argument("-t", default='1', metavar='', help=('''
number of CPUs
'''
))
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

# assign parameters
cdhit_file = args.cdhit
copy_file = args.fasta
out_dir = args.o
# max_size = args.max_size


if out_dir[-1] == "/":
    out_dir = out_dir[:-1]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
out_dir = os.path.abspath(out_dir)

n_cpus = int(args.t)


def parse_cdhit(cdhit):
    """
    Parse cd-hit-est output file having clusters. It returns a dict where the key is the representative sequence ID of the
    cluster and the value the list of other sequences ID in the cluster.
    """

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
    """
    Given a fasta file, it returns a dict where the key is the ID of the sequence and the value is the sequence itself.
    It can be used with the fasta file used for cd-hit-est.
    """

    seqs = {} # key: name of the copy, value: sequence
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name = seq_record.id
        seq = seq_record.seq
        seqs[name] = str(seq).upper() # it allows to avoid errors while running Refiner
    return seqs


def refiner(cluster_file):
    """
    Call Refiner in a subprocess.
    """
    command = ['Refiner', cluster_file] 
    subprocess.run(command)


def do_work(ncpu, out_dir):
    """
    Recovers all fasta files for cluster present at the moment in the folder and call Refiner command using parallelization.
    """

    tmp_clst = glob.glob(out_dir + '/*.clst.fasta') # select all fasta files for each cluster present in the directory
    
    pool = mp.Pool(ncpu)
    for clst in tmp_clst:
        pool.apply_async(refiner, (clst, ))              

    pool.close()
    pool.join()

    for clst in tmp_clst:
        os.remove(clst)


def update_and_clean(out_dir):
    """
    Add new consensi to consensi.fa and consensi.stk and clean temporary files.
    """
    # output files
    consensi = out_dir + "/consensi.fa"
    stk = out_dir + "/consensi.stk"

    # add the sequence to consensi.fa
    tmp_cons = glob.glob(out_dir + '/*.refiner_cons')

    with open(consensi, "a+") as f:
        for cons_fasta in tmp_cons:
            with open(cons_fasta, "r") as ref_file:
                data = ref_file.read()
                f.write(data)
            os.remove(cons_fasta)
                
    # add the alignement to consensi.stk
    tmp_stk = glob.glob(out_dir + '/*.refiner.stk')

    with open(stk, "a+") as f:
        for cons_stk in tmp_stk:
            with open(cons_stk, "r") as ref_file:
                data = ref_file.read()
                f.write(data)
            os.remove(cons_stk)

    tmp_RM = glob.glob(out_dir + '/RM_*')
    for rm in tmp_RM:
        path = os.path.join(out_dir, rm)
        if os.path.isdir(path):
            shutil.rmtree(path)


def plot_dist(cdhit_dict, out_dir):
    """
    Plot the distribution (with an histogram) of the number of sequences in clusters.
    """

    dist = []
    for k in cdhit_dict:
        dist.append(len(cdhit_dict[k]) + 1)
    # x="Number of sequences in the cluster",
    fig = px.histogram(dist, marginal="box", title="Distribution of the number of sequences")
    fig.write_html(out_dir + "/dist.html")
    

def call_consensus(cdhit_dict, seqs_dict, out_dir, ncpu):
    """
    It separates singleton clusters in singleton.fa and then it call a consensus for each cluster using Refiner (from 
    RepeatModeler2) in consensi.fa and the alignement at .stk format in consensi.stk.
    It uses dict generated by parse_cdhit and parse_fasta functions. 
    out_dir: output directory.
    ncpu: number of precessor to use for multiprocessing.
    """
    
    # counters for statistics
    nb_repr_within_copies = 0
    nb_singleton = 0
    nb_low_copy = 0
    nb_cluster = 0

    # output files
    singleton = out_dir + "/singleton.fa"
    low_copy = out_dir + "/low_copy_number.fa" # copy number < 5

    c = 0 # counter for naming consensus sequences

    # main   
    nb_cluster_file = 0 # to count how many file generate for parallelization
    os.chdir(out_dir)

    for seq_id in seqs_dict.keys():

        if seq_id in cdhit_dict.keys():
            nb_repr_within_copies += 1

            if len(cdhit_dict[seq_id]) == 0: # the cluster is singleton
                header = str(seq_id).replace("|", "_")
                nb_singleton += 1                
                with open(singleton, "a+") as f:                    
                    f.write(">" + header + "\n" + str(seqs_dict[seq_id]) + "\n")

            elif len(cdhit_dict[seq_id]) <= 5: # the cluster is a low copy family
                
                seqs_id = cdhit_dict[seq_id]
                seqs_id.append(seq_id) # add representative sequence

                nb_low_copy += 1
                
                with open(low_copy, "a+") as f:
                    for seq in seqs_id:
                        f.write(">" + str(seq) + "\n" + seqs_dict[seq] + "\n")
                    f.write("//\n") # to separate different clusters


            else: # the cluster is not singleton or low copy
                
                nb_cluster += 1
                nb_cluster_file +=  1

                cluster_file = out_dir + "/" + str(nb_cluster_file) + ".clst.fasta"
                seqs_id = cdhit_dict[seq_id]
                seqs_id.append(seq_id) # add representative sequence
                
                # write sequences of the cluster in a fasta file
                with open(cluster_file, "w") as f:
                    c = 0
                    for seq in seqs_id:
                        f.write(">" + "sequence" + str(c) + "\n" + seqs_dict[seq] + "\n")
                        c += 1
                        if c == 500:
                            break
                                
                if nb_cluster_file == ncpu:
                    do_work(ncpu, out_dir)
                    update_and_clean(out_dir)    
                    nb_cluster_file = 0                        

                c += 1

    if nb_cluster_file != 0:
        do_work(nb_cluster_file)
        update_and_clean(out_dir)  

    with open(out_dir+"/stats.txt", "w") as stats:
        stats.write(f"""
        Total number of copies:\t{str(len(seqs_dict))}
        Number of cluster:\t{str(len(cdhit_dict))}
        Number of representative sequences founded in the fasta file with copies:\t{str(nb_repr_within_copies)}
        Number of singleton (clusters with only one sequence):\t{str(nb_singleton)}
        Number of cluster with a low number of copies (<= 5):\t{str(nb_low_copy)}
        Number of cluster with 6 or more sequences:\t{str(nb_cluster)}
        """)
    
    # plot_dist(cdhit_dict, out_dir)



if __name__ == "__main__":
    start = time.time()
    
    clusters = parse_cdhit(cdhit_file)
    seqs = parse_fasta(copy_file)
    call_consensus(clusters, seqs, out_dir, n_cpus)
    
    end = time.time()
    elapsed = end - start
    print(f"\n\n\nElapsed time: {elapsed}")