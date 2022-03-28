#!/bin/bash

## Author: Tommaso Barberis

## Date: 22/03/2022

## Description: script that allows to parse OneCodeToFindThemAll output and recover a fasta file containing sequences from of all transposable elements copies

## Parameters:
# -i|--input: directory that contain output from OneCodeToFindThemAll; OPTIONAL, it can detect automatically the output from OneCodeToFindThemAll if the script is executed in the same folder
# -g|--genome: path the reference genome used to call RepeatMasker; REQUIRED

## Dependencies:
# python3 (tested with python v.3.7.9) and the following modules:
#	argparse
#	os
# bedtools (tested with v.2.30.0)


start_time=$(date +%s)
set -e

POSITIONAL_ARGS=()

INPUT_DIR=./
SCRIPTS=~/scripts

## Parsing parameter
while [[ $# -gt 0 ]]; do
	case $1 in
		-i|--input)
			INPUT_DIR="$2"
			shift
			shift
			;;
		-g|--genome)
			GENOME="$2"
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

## input check
#INPUT_DIR=$(realpath $INPUT_DIR 2> /dev/null)
[[ ! -d $(realpath $INPUT_DIR) ]] && { echo -e "Invalid input directory!"; exit 1; }
GENOME=$(realpath $GENOME 2> /dev/null) || { echo -e "Genome file is not provided!"; exit 1; }
[[ ! -f $GENOME ]] && { echo -e "Invalid genome file!"; exit 1; } #TODO: check if it is a FASTA file

## dependencies check
dependencies=(python3 bedtools) #TODO: check if it works
for dep in ${dependencies[@]}; do
	if command -v $dep > /dev/null 2>&1; then
		true
	else
		echo -e "$dep not found!"
		exit 1
	fi
done

## main function

echo -e "## Summarising octfta results in a single .tsv file (recovering only full-length copies)"

# recover elements from each contig/scaffold/chromosome in final.tsv
cat *elem_sorted.csv | grep -v "Score" > final.tsv

# recover only full copies
HEAD=$(echo -e "$(head -1 $(ls -1 *sorted* | head -1))")
echo $HEAD > final.full.tsv
grep "###" final.tsv >> final.full.tsv

# convert to bed file
echo -e "## Converting the .tsv file to bed file"
python3 $SCRIPTS/octfta_to_bed.py -octfta final.full.tsv

# recover sequences from full copies
echo -e "## Recovering sequences of copies in a FASTA file"
bedtools getfasta -fi $GENOME -bed final.full.bed -name > copies.fasta

## Computing time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
#format=$(date -ud "@$elapsed" +"$(expr s% / 3600 / 24 ) days %H hr %M min %S sec")
#echo -e "Elapsed time: $format"
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
