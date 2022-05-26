"""
Author: Tommaso Barberis
Date: 03/05/2022
Description: artificial trimming of the ends of sequences
"""

from Bio import SeqIO
import argparse, os

def trim_ends(fasta_file, thres):
    
    fields = fasta_file.split(".")
    prefix = ""
    for x in fields[:-1]:
        prefix = prefix + "." + str(x)
    prefix = prefix[1:]
    
    new_name = prefix + ".trim." + fields[1]

    if os.path.exists(new_name):
        print(f"A trimmed file ({new_name}) already exists")
        return None

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        
        seq_len = len(str(seq_record.seq))
        trim_len = int(seq_len * thres)
        new_seq = str(seq_record.seq)[trim_len:-trim_len]
        
        with open(new_name, "a+") as f:
            f.write(">" + seq_record.name + "\n" + new_seq + "\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", metavar='', help=(f'''
    FASTA file
    '''))
    parser.add_argument("-t", metavar='', default='0.1', help=(f'''
    threshold to cut at end, default it cuts 10% (0.1) of the length on each end.
    '''))

    args = parser.parse_args()

    # assign parameters
    fasta_file = args.f
    thres = float(args.t)

    trim_ends(fasta_file, thres)