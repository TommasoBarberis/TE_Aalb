# -*-coding:Utf-8 -*

"""
Author: Tommaso Barberis
Date: 22/03/2022
Description: convert tab-delimited file containing full elements recovered using OneCodeToFindThemAll in a .bed file
"""

import argparse


parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("-octfta", help=('''
.tsv input file from OneCodeToFindThemAll containing only full copies (lines that start with "###").
 '''))
group.add_argument("-bed", help=('''
.bed file to convert as OCTFTA output but with a supplementar column on the left for fragments.     
 '''))


args = parser.parse_args()

octfta_file = args.octfta
bed_file = args.bed

def octfta_to_bed(octfta):
    output_file = octfta.replace(".tsv", "") + ".bed"
    with open(octfta, "r") as f:
        header = f.readline() # skip header

        lines = f.readlines()
        for line in lines:
            line = line.split()

            chrom = line[4] # Chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671) name (REQUIRED)
            chrom_start = line[5] # Start coordinate on the chromosome or scaffold for the sequence considered (the first base on the chromosome is numbered 0) (REQUIRED)
            chrom_end = line[6] # End coordinate on the chromosome or scaffold for the sequence considered. This position is non-inclusive, unlike chromStart. (REQUIRED)
            name = line[9] # Name of the line in the BED file, name of the element in the OneCodeToFindThemAll
            score = line[0] # Score between 0 and 1000
            strand = line[8] # DNA strand orientation (positive ["+"] or negative ["-"] or "." if no strand)
            subs_pourcentage = line[1]
            del_pourcentage = line[2]
            ins_pourcentage = line[3]
            family = line[10]
            rm_ID = line[14]
            number_hits = line[15]

            with open(output_file, "a") as o:
                new_line = chrom + "\t" + chrom_start + "\t" + chrom_end + "\t" + name + "\t" + score + "\t" + strand + \
                    "\t" + subs_pourcentage + "\t" + del_pourcentage + "\t" + ins_pourcentage + "\t" + family + "\t" + \
                    rm_ID + "\t" + number_hits + "\n"
                o.write(new_line)
                

def bed_to_octfta():
    pass
        
        
if __name__ == "__main__":
    if octfta_file is not None:
        octfta_to_bed(octfta_file)
    if bed_file is not None:
        bed_to_octfta(bed_file)