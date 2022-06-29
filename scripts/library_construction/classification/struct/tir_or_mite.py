"""
Author: Tommaso BARBERIS
Date: 26/05/2022
Description: split putative tir in TIR or in MITE by the length
"""

import argparse
from Bio import SeqIO


def tir_or_mite(fasta_file, seq_file):
    with open(seq_file, 'r') as f:
        seq_ids = [x.replace("\n", "") for x in f.readlines()]

    tir = {}
    mite = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = seq_record.id
        if seq_id in seq_ids:
            seq = str(seq_record.seq)
            if len(seq) <= 750:
                mite[seq_id] = seq
            else:
                tir[seq_id] = seq
    
    with open("mite.fa", "w") as f:
        for s in mite.keys():
            f.write(">" + s + "\n" + mite[s] + "\n")
    with open("tir.fa", "w") as f:
        for s in tir.keys():
            f.write(">" + s + "\n" + tir[s] + "\n")
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", metavar='', help=('''
    FASTA file with putative sequences
    '''))
    parser.add_argument("-t", metavar='', help=('''
    text file with sequences ID identified by einverted
    '''))

    args = parser.parse_args()

    # assign parameters
    fasta_file = args.f
    seq_file = args.t

    if fasta_file and seq_file:
        tir_or_mite(fasta_file, seq_file)
    pass