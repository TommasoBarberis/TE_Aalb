#!/bin/bash

# Author: Tommaso Barberis
# Date: 22/04/2022
# Description: test for end extension of consensi

# $1: fasta file with consensi 
# $2: ref genome
# $3: cpus

# dependancies:
# -mafft
# -scripts from https://github.com/annaprotasio/TE_ManAnnot/blob/main/bin

source ~/anaconda3/etc/profile.d/conda.sh
conda activate TE_Aalb
export fasta_file=$(realpath $1)

work_dir=$(realpath .)
genome_file=$(realpath $2)

python3 - << EOF
from Bio import SeqIO
import os

c = 0
for seq_record in SeqIO.parse("$fasta_file", "fasta"):
    os.mkdir("seq_" + str(c))

    with open("seq_" + str(c) + "/seq.fasta", "w") as f:
        f.write(">sequence\n" + str(seq_record.seq))
    c += 1
EOF

for dir in $(ls -d -1 $work_dir/seq_*); do

    echo -e "\npreparatory stuff \n $dir \n\n"
    cd $(realpath $dir)
    make_fasta_from_blast.sh $genome_file seq.fasta 80 1500
    nb_seq=$(grep -c ">" seq.fasta.blast.bed.fa)
    if [ $nb_seq -gt 200 ]; then
        ready_for_MSA.sh seq.fasta.blast.bed.fa 200 50
        mafft --reorder --thread $3 seq.fasta.blast.bed.fa.rdmSubset.fa > maf.fa
    else
        mafft --reorder --thread $3 seq.fasta.blast.bed.fa > maf.fa
    fi
    
    echo -e "\n CIAlign stuff \n\n"
    CIAlign --infile maf.fa --outfile_stem no_clean --make_consensus --plot_output --plot_format svg
    CIAlign --infile maf.fa --outfile_stem no_ins --remove_insertions --make_consensus --plot_output --plot_format svg

done

conda activate te_aid

for dir in $(ls -d -1 $work_dir/seq_*); do
    cd $(realpath $dir)
    TE-Aid -g $genome_file -q seq.fasta
    TE-Aid -g $genome_file -q no_clean_consensus.fasta
    TE-Aid -g $genome_file -q no_ins_consensus.fasta
done