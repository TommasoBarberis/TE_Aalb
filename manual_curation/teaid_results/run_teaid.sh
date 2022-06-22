#!/bin/bash

## Description: run TE-Aid on fasta sequences contained in a file

# $1: fasta file
# $2: genome file


while read -r header; do
    read -r sequence
    

    file_name=$(echo $header | sed 's/>//g')
    
    # create tmp file with only one sequence
    echo -e "$header\n$sequence" > $file_name.fa
    TE-Aid -q $file_name.fa -g $2
    rm $file_name.fa
done < $1