"""
Author: Tommaso BARBERIS
Date: 05/05/2022
Description: python wrapper for polishTE (https://github.com/TommasoBarberis/TE_Aalb/blob/main/scripts/library_construction/extender) allowing multiprocessing
             and  to processing multifasta
"""

import argparse, os, shutil, subprocess
from Bio import SeqIO
import multiprocessing as mp
import plotly.express as px
    

def run_polish(c, seq_record, results_file, genome):
    """
    Run polishTE
    """
    
    seq_dir = "polishte_tmp/seq_" + str(c)
    os.mkdir(seq_dir)
    
    seq_filename = "polishte_tmp/seq_" + str(c) + "/seq.fasta"
    
    with open(seq_filename, "w") as f:
        f.write(">" + seq_record.description + "\n" + str(seq_record.seq) + "\n")
    
    polish_cmd = ['extender', '-i', seq_filename, '-g', genome, '-t', '1', '-ins', '200', '-o', seq_dir]
    subprocess.run(polish_cmd)

    def get_seq(new_file):
        for seq_new_file in SeqIO.parse(new_file, "fasta"):
            seq = str(seq_new_file.seq)
        
        return seq

    if os.path.exists(seq_dir + "/new_consensus.fasta"):
        seq = get_seq(seq_dir + "/new_consensus.fasta")
    elif os.path.exists(seq_dir + "/simple_repeat.fasta"):
        seq = get_seq(seq_dir + "/simple_repeat.fasta")
    else:
        seq = str(seq_record.seq)

    with open(results_file, "a+") as r:
        r.write(">" + seq_record.description + "\n" + seq + "\n")

    # for stats
    if len(seq) > len(str(seq_record.seq)):
        diff = len(seq) - len(str(seq_record.seq))
        stat = "extend"
    elif len(seq) < len(str(seq_record.seq)):
        diff = len(str(seq_record.seq)) - len(seq)
        stat = "trim"
    else:
        diff = 0
        stat = "NA"

    percent = (100 * diff) / len(str(seq_record.seq))
    with open("stat.txt", "a+") as s:
        s.write(stat + "\t" + str(percent) + "\n")

    shutil.rmtree(seq_dir)

def main(fasta_file, genome, ncpu):
    """
    Manage FASTA sequences
    """

    # if a results file already exists, quit the program
    fields = fasta_file.split(".")
    results_file = ".".join(fields[:-1])
    results_file = ".".join([results_file, "polished", fields[-1]])
    if os.path.exists(results_file):
        print("A results file already exists")
        quit()

    # init the multiprocessing pool
    pool = mp.Pool(processes=ncpu)

    # crate tmp folder
    try:
        os.mkdir("polishte_tmp")
    except:
        shutil.rmtree("polishte_tmp")
        os.mkdir("polishte_tmp")

    # run polishTE on each sequence
    c = 1
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        pool.apply_async(run_polish, (c, seq_record, results_file, genome))
        
        c += 1

        if c == 5:
            break

    pool.close()    
    pool.join()

    # plot stats
    extend_percent_lst = []
    trim_percent_lst = []
    with open("stat.txt", "r") as s:
        data = s.readlines()
        for line in data:
            fields = line.split()
            if fields[0] == "extend":
                extend_percent_lst.append(float(ields[1]))
            elif fields[0] == "trim":
                trim_percent_lst.append(float(fields[1]))


    fig = px.histogram(extend_percent_lst, marginal="box", title="Extension percent")
    fig.update_layout(
        xaxis_title="percent", 
        showlegend=False,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    fig.write_image("extend.pdf")

    fig1 = px.histogram(trim_percent_lst, marginal="box", title="Trim percent")
    fig1.update_layout(
        xaxis_title="percent", 
        showlegend=False,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )    
    fig1.write_image("trim.pdf")

    # clean
    os.remove("stat.txt")
    shutil.rmtree("polishte_tmp")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", metavar='', help=('''
    FASTA file with consensi to polish
    '''))
    parser.add_argument("-g", metavar='', help=('''
    genome file
    '''))
    parser.add_argument("-t", metavar='', default='1', help=('''
    number of cpus
    '''))

    args = parser.parse_args()

    # assign parameters
    fasta_file = args.f
    genome = args.g
    cpus = int(args.t)

    main(fasta_file, genome, cpus)