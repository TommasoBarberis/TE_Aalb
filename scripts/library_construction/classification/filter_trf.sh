#!/bin/bash

# Author: Tommaso Barberis
# Date: 20/06/2022
# Description: script to filter trf results (based on TRFMask script from RepeatModeler2)
# Parameters:
#   $1 -> results from trf (runned with: trf4 file.fa 2 7 7 80 10 50 500 -ngs -h > tandem_repeats.txt)


while read -r line; do
    
    if [[ $line == @* ]]; then
        if [ -f repeats ]; then
            
            echo $id >> $1.filtered
            cat repeats >> $1.filtered
            rm repeats
        fi
        
        id=$line

    else
        period=$(echo $line | awk '{print $3}')
        copies=$(echo $line | awk '{print $4}')
        entropy=$(echo $line | awk '{print $13}')

        if (( $(echo "$period > 1" | bc -l) )) && (( $(echo "$copies > 8" | bc -l) )) || (( $(echo "$entropy > 1.95" | bc -l) )); then
                echo $line >> repeats
        fi
    fi
done < $1