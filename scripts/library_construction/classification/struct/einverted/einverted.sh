#!/bin/bash

# Author: Tommaso BARBERIS
# Date: 20/06/2022
# Description: script to run einverted
# Parameters:
#   $1: FASTA file

while read -r header; do
    read -r seq

    id=$(echo $header | sed 's/>//g')
    echo $id
    echo -e "$header\n$seq" > $id.fa
    
    einverted $id.fa -gap 12 -threshold 50 -match 3 \
        -mismatch -4 -outfile $id.inv -outseq $id.out.fasta \
        -sformat1 fasta

    cat $id.inv >> putative.inv
    cat $id.out.fasta >> putative.out.fasta

    rm $id.fa $id.inv $id.out.fasta

done < $1
