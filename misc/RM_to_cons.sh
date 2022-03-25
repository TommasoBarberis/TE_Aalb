#!/bin/bash

## Authors:
# Daniel Siqueira de Oliveira => v.1 feb/2022 clustering_cons.sh
# Tommaso Barberis => v.2 feb/2022 RM_to_cons.sh

## Description:
# Parse .gff file obtained with -gff option in RepeatMasker and return consensus sequences

## Dependencies: all dependencies have to stay in the global PATH
# tee 
# bedtools getfasta (tested with 2.30.0)
# cd-hit-est (tested with v.4.8.1)
# conda env: TE_Aalb
# Refiner (from RepeatModeler2, tested with v2.0.2)

## Parameters:
# -a|--gff: gff file from repeatMasker output; REQUIRED
# -g|--reference: genome reference at the fasta format used to run RepatMasker; REQUIRED
# -o|--output: output dir; OPTIONAL, default=./
# --overwrite: ecrase old job in the same directory; OPTIONAL, default=NO
# -t|--cpu: number of threads to use; OPTIONAL, default=1

start_time=$(date +%s)

POSITIONAL_ARGS=()                      

OUTPUT=./
OVERWRITE=NO
CPU=1

## Parsing parameter
while [[ $# -gt 0 ]]; do
	case $1 in
		-a|--gff)
			GFF="$2"
			shift
			shift
			;;
		-g|--reference)
			REF="$2"
			shift
			shift
			;;
		-o|--output)
			OUTPUT="$2"
			shift
			shift
			;;
		-t|--cpu)
			CPU="$2"
			shift
			shift
			;;
		--overwrite)
			OVERWRITE=YES
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

## Log file
log_time=$(date '+%d-%m-%Y_%H:%M:%S')
log_file=$OUTPUT/RMtoCons-$log_time.log

log_msg () {
	if [ "$2" = "only_log" ];then
		echo -e "$1" >> $log_file
	else
		echo -e "$1" 2>&1 | tee -a $log_file
	fi
}

## Error handling

exit_msg () {
	echo -e "line $LINENO: ERROR! $1" 2>&1 | tee -a $log_file
	 exit 1
}

## Input check

log_msg "## Input check"
gff_file=$(realpath $GFF 2> /dev/null) || exit_msg ".gff file not found!"
genome_file=$(realpath $REF 2> /dev/null) || exit_msg "genome file not found!"
output_dir=$(realpath $OUTPUT 2> /dev/null) || exit_msg "output directory not valid!"
if [ -n "$CPU" ] && [ "$CPU" -eq "$CPU" ] 2>/dev/null; then
	true
else
	exit_msg "invalid CPUs option!"
fi

## Check dependencies installation

log_msg "## Dependencies check"
dependencies=(tee bedtools cd-hit-est Refiner)
for dep in ${dependencies[@]}; do
	if command -v $dep >/dev/null 2>&1 ; then
		log_msg "\t $dep found" only_log
	else
		exit_msg "$dep not found"
	fi
done

conda_envs=$(conda env list | grep -v "#" | grep -v "^$" | awk '{print $1}') || exit_msg "conda not installed"
if { echo $conda_envs | grep 'TE_Aalb'; } >/dev/null 2>&1; then
	log_msg "TE_Aalb env found" only_log
else
	exit_msg "TE_Aalb env not installed"
fi

#TODO: add OCTFTA part
## Parsing TE copies from genome reference

if [ ! -f TE_copies.fasta ]; then	
	log_msg "## Parsing TE copies from genome reference" 2>&1
	bedtools getfasta -fi $genome_file -bed $gff_file > TE_copies.fasta
else
	log_msg "## Parsing TE copies from genome reference, SKIP: copies.fasta already present"
fi
seq_fasta=$(realpath TE_copies.fasta)

## Clustering sequences

if [ ! -f TE_copies.fasta.hit.clstr ]; then
	log_msg "## Clustering sequences" 
	cd-hit-est -i TE_copies.fasta -o TE_copies.fasta.hit -c 0.8 -n 4 -M 10000 -T $CPU -d 0
else
	log_msg "## Clustering sequences, SKIP: TE_copies.fasta.hit.clstr already present"
fi
clust=$(realpath TE_copies.fasta.hit.clstr)

## Call consensi

log_msg "## Call consensi from each cluster"
conda activate TE_Aalb
python3 ~/scripts/parse_cdhit.py -i $clust -fasta $seq_fasta -t $CPU

## Computing time
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
format=$(date -ud "@$elapsed" +"$(expr %s / 3600 / 24 ) days %H hr %M min %S sec")
log_msg "\n\n Elapsed time: $format"
