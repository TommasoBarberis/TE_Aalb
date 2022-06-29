#!/bin/bash

# Author: Tommaso BARBERIS
# Date: 27/05/2022
# Description: parse fasta sequence with classification

# $1 fasta file
# $2 file with classification

while read -r line; do
    seq_id=`echo $line | awk '{print $1}'`
    sf=`echo $line | awk '{print $2}'`
    fam=`echo $line | awk '{print $3}'`

    if [ $fam == '-' ]; then
        header=$(echo ">$seq_id#$sf")
    else
        header=$(echo ">$seq_id#$sf/$fam")
    fi

    seq=`grep -A1 "$seq_id$" $1 | tail -1`
    
    echo -e "$header\n$seq" >> $1.classified

done < $2