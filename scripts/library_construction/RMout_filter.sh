#!/bin/bash

## Author: Tommaso Barberis 

## Date: 25/03/2022

## Description: custom script to filter RepeatMasker output.
# It allows to filter hits by length (<80bp).

## Dependancies:


## Parameters:
#-i|--input: RepeatMasker .out file; REQUIRED
# -o|--output: output dir (default: ./); OPTIONAL
# -l|--min-length: minimal length in bp for the hit, based on position in the query (default: 80), if -l 0, it skips this step; OPTIONAL
# -a|--annotation: file with key words to filter out, based on columns 'matching repeats' and 'repeat class/family'; OPTIONAL

POSITIONAL_ARGS=()

OUTPUT=./
LENGTH=80

## Parsing parameter
while [[ $# -gt 0 ]]; do
	case $1 in
		-i|--input)
			INPUT="$2"
			shift
			shift
			;;
		-o|--output)
			OUTPUT="$2"
			shift
			shift
			;;
        -l|--min-length)
			LENGTH="$2"
			shift
			shift
			;;
        -a|--annotation)
            ANN="$2"
            shift
            shift
            ;;
		-*|--*)
			echo "Unknown option $1"
			exit 1
			;;
		*)
			POSITIONAL_ARGS+=("$1")
			shift
			;;
	esac
done

rm_file=$(realpath $INPUT)

OUTPUT=$(realpath $OUTPUT)

filtered_file=${rm_file##*/}
filtered_file=$(echo "$OUTPUT/$filtered_file.filtered")
rm $filtered_file 2> /dev/null


# recover the header
head -3 $rm_file >> $filtered_file

# filter by length
if [ $LENGTH = 0 ]; then
    true
else
    tail -n +3 $rm_file | while read -r line; do
        awk -v len=$LENGTH '{x=$7-$6; if (x>len) print $0}' >> tmp
    done 
fi

# filter by annotation
if [  "$ANN" = '' ]; then
    true # skip this step
else
    ann_file=$(realpath $ANN)
    
    if [ -f $ann_file ]; then
        if [ -f tmp ]; then
            grep -v -f $ann_file tmp >> tmp1
        else
            tail -n +3 $rm_file | grep -v -f $ann_file >> tmp1
        fi
    else
        echo "Annotation file not found!"
        exit 1
    fi    
fi

# write results
if [ -f tmp1 ]; then
    cat tmp1 >> $filtered_file
elif [ -f tmp ] && [ ! -f tmp1 ]; then
    cat tmp >> $filtered_file
fi

rm tmp tmp1 2> /dev/null