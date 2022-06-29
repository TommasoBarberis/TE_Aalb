#!/bin/bash

# Author: Tommaso Barberis
# Date: 20/06/2022
# Description: script to filter trf results (based on TRFMask script from RepeatModeler2)
# Parameters:
#   $1 -> results from trf (runned with: trf4 file.fa 2 7 7 80 10 50 500 -ngs -h > tandem_repeats.txt)
#   $2 -> minum number of copies (suggested: 5)
#   $3 -> FASTA file (file used with TRF)
#   $4 -> threshold on sequence coverage (suggested: 0.8)


while read -r line; do
    
    if [[ $line == @* ]]; then
        if [ -f repeats ]; then
            
            if (( $(echo "$cov > $4*$seq_len" | bc -l) )); then
                echo $id >> $1.filtered
                cat repeats >> $1.filtered
            fi
            rm repeats
        fi
        
        id=$line
        seq=$(grep -A1 "$id$" $3 | tail -1)
        seq_len=$(echo $seq | wc -c)
        cov=0

    else
        start=$(echo $line | awk '{print $1}')
        end=$(echo $line | awk '{print $2}')
        copies=$(echo $line | awk '{print $4}')

        if (( $(echo "$copies > $2" | bc -l) )); then
        #  && (( $(echo "$end - $start > $4*$seq_len" | bc -l) ))
                echo $line >> repeats
                rep_cov=$(($end-$start))
                cov=$(($cov+$rep_cov))
        fi
    fi
done < $1