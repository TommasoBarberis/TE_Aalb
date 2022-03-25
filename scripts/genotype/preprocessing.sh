#!/bin/bash

## Author: Tommaso Barberis (inspired from https://imbs-hl.github.io/illumina_seq.html)

## Description: script to pre-processing raw reads for TE insertion calling

## Dependancies:
# tee: for logging
# fastqc: quality control check, tested on FastQC - 0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip) 

## Parameters:
# -f|--file = fastq file; REQUIRED
# -o|--output = output directory (default: ./); OPTIONAL
# --overwrite = overwrite older jobs; OPTIONAL
# -phix|--phix-genome = PATH to PhiX genome; REQUIRED
# genome reference dir (indexed with bwa-mem2)

POSITIONAL_ARGS=()

OUTPUT=./
OVERWRITE=NO

## Parsing parameter
while [[ $# -gt 0 ]]; do
	case $1 in
		-f|--file)
			FILE="$2"
			shift
			shift
			;;
		-o|--output)
			OUTPUT="$2"
			shift
			shift
			;;
		--overwrite)
			OVERWRITE=YES
			shift
			;;
		-phix|--phix-genome)
			PHIX_GENOME="$2"
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

## Log file

log_time=$(date '+%d-%m-%Y_%H:%M:%S')
log_file=$OUTPUT/preprocessing-$log_time.log


## Error handling

exit_msg () {
	echo "line $LINENO: ERROR! $1"
	exit 1
}


## input check

filename=${FILE##*/}
filename=${filename%.*}

FILE=$(realpath $FILE)
OUTPUT=$(realpath $OUTPUT)
reference_file_phix=$(realpath $PHIX_GENOME 2> /dev/null) || exit_msg "PhiX genome not found!"

## create directory to pre-process the sample

work_dir=$OUTPUT/pre-processing/$filename

if [ ! -d $work_dir ] || [ $OVERWRITE == "YES" ]; then
	[[ -d $work_dir ]] && rm -r $work_dir
	mkdir -p $OUTPUT/pre-processing/$filename
	echo -e "## Creation of the output folder for pre-processing" 2>&1 | tee -a $log_file

	cd $work_dir

	## First quality control
	if [ ! -d $work_dir/QC1/ ]; then
		mkdir -p $work_dir/QC1/
		echo -e "## Quality control on raw reads" 2>&1 | tee -a $log_file
		fastqc $FILE -o $work_dir/QC1
	fi

	## Detection of contamination
	if [ ! -d $work_dir/contam/ ]; then
		echo -e "## Detection of contaminations" 2>&1 | tee -a $log_file
		echo $reference_file_phix
	fi

else
	echo -e "## Pre-processing of $filename already done" 2>&1 | tee -a $log_file
fi

gatk_resources_dir=/beegfs/data/tbarberis/TE_alb/genotyping/pre-processing/contam/gatk_bundle
known_sites_file_gold=$gatk_resources_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf
known_sites_file_1000g=$gatk_resources_dir/1000G_phase1.snps.high_confidence.hg38.vcf
known_sites_file_dbsnp=$gatk_resources_dir/Homo_sapiens_assembly38.dbsnp138.vcf
known_sites_file_hapmap=$gatk_resources_dir/hapmap_3.3.hg38.vcf
known_sites_file_omni=$gatk_resources_dir/1000G_omni2.5.hg38.vcf

#rg_id=<FID>.<Lx>
#rg_sm=<SID>
#rg_lb=<LIB>
#rg_pu=<FID>.<Lx>.<SID>
#rg_pl="ILLUMINA"
#read_group="@RG\tID:$rg_id\tSM:$rg_sm\tLB:$rg_lb\tPU:$rg_pu\tPL:$rg_pl"
