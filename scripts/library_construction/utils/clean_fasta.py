"""
Author: Tommaso Barberis
Date: 12/04/2022
Description: clean FASTA file from sequences with a given character
"""

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("-i", metavar='', help=('''
input FASTA file
'''))
parser.add_argument("-o", metavar='', help=('''
output FASTA file
'''))
parser.add_argument("-c", metavar='', help=('''
character to purge (case sensitive)
'''))

args = parser.parse_args()

input_file = args.i 
output_file = args.o
char = args.c

def main(input_f, output_f, char):
    for seq_record in SeqIO.parse(input_f, "fasta"):
        seq_name = str(seq_record.description)
        seq = str(seq_record.seq)
        if char in seq:
            pass
        else:
            with open(output_f, "a+") as f:
                f.write(">" + seq_name + "\n" + seq + "\n")
    pass

if __name__ == "__main__":
    main(input_file, output_file, char)