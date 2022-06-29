#!/bin/bash

# Author: Tommaso BARBERIS
# Date: 20/06/2022
# Description: script to filter einverted output (FASTA file) in order 
#   to obtain sequences with TIR in first and last n nucleotides
# Parameters:
#   $1: FASTA file (output from einverted)
#   $2: threshold, nb of nucleotide from start and from end to contain the TIR (suggested: 100nt) 
#   $3: FASTA file with original sequences (file used with einverted)
#   $4: threshold to separate TIR elements from MITE elements (suggested: 800nt)


grep ">" $1 > filter_tmp

while read -r tir_5; do
    read -r tir_3

    id=$(echo $tir_5 | sed 's/>//g' | awk -F"_" '{print $1"_"$2}')
    seq=$(grep -A1 "$id$" $3 | tail -1)
    seq_len=$(echo $seq | wc -c)
    
    tir_5_end=$(echo $tir_5 | sed 's/>//g' | cut -f4 -d"_")
    tir_3_begin=$(echo $tir_3 | sed 's/>//g' | cut -f3 -d"_")
    
    if [[ $tir_5_end -le $tir_3_begin ]] && [[ $tir_5_end -le $2 ]] && (( $(echo "$tir_3_begin > ($seq_len - $2)" | bc -l ) )); then
        if (( $(echo "$seq_len > $4" | bc -l ) )); then
            echo $id >> tir.lst
        else
            if (( $(echo "$seq_len > $2*2" | bc -l) )); then
                echo $id >> mite.lst            
            fi
        fi
    fi


done < filter_tmp

rm filter_tmp
