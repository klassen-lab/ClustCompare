#!/usr/bin/bash

# Script to annotate a batch of NCBI genomes using antiSMASH. Assumes that run_antismash.py works from terminal
# Will put annotations in ../Data/ and modify ../CHANGELOG.txt
# v3.0 January 2, 2018 - adapted to Jan 2018 script versions. Requires all genomes be located in ../Data/genomes. Runs on .gbk, .gbff, .fasta, .fna, .fa nucleotide formats.
# v3.1 February 15, 2018 - Accommodates genomes.list output file
# v3.2 February 19, 2018 - adds command line options for number of threads and data paths
# v3.3 February 21, 2018 - modified data paths, accommodated new cluster_renamer.pl

usage()
{
	echo "USAGE:"
	echo "Syntax: ./antiSMASH_annotation.sh [OPTIONS]"
	echo "Options:"
	echo "    -c    number of CPU threads"
	echo "    -i    path to directory containing input data files"
	echo "    -o    path to directory where output data files should be placed"
	echo "    -r    path to directory where output results files should be placed"
	echo "    -l    path to log file to update"
	echo "    -h    print this help message"
}

# default parameters
NPROC=$(nproc)
INPUT_DATA_PATH=../Data/genomes
OUTPUT_DIR=../Data/antiSMASH_annotations
OUTPUT_RESULTS_DIR=../Results
LOG_FILE=../CHANGELOG.txt

# parse command line options
while [ "$1" != "" ]; do
	case $1 in
		-c )	shift
			NPROC=$1
			;;
		-i )	shift
			INPUT_DATA_PATH=$1
			;;
		-o )	shift
			OUTPUT_DIR=$1
			;;
		-r )	shift
			OUTPUT_RESULTS_DIR=$1
			;;
		-l )	shift
			LOG_FILE=$1
			;;
		-h )	usage
			exit
			;;
		* )	usage
			exit 1
	esac
	shift
done

echo "Using $NPROC cores"
echo "Looking for input data here: $INPUT_DATA_PATH"
echo "Ouput data files will be here: $OUTPUT_DIR"
echo "Ouput data results files will be here: $OUTPUT_RESULTS_DIR"
echo "Updating log: $LOG_FILE"

# Add note to log
echo -e "## `date +%Y-%m-%d:%H:%M:%S` \n\n* antiSMASH_annotations.sh v3.3 started" >> $LOG_FILE 


# check for $INPUT_DATA_PATH directory
if [ -d $INPUT_DATA_PATH ]; then
	echo "$INPUT_DATA_PATH found"
else
	echo "$INPUT_DATA_PATH not found - exiting antiSMASH_annotation.sh"
	exit
fi


# unzip any zipped files, if they exist
for f in $(eval echo "$INPUT_DATA_PATH/*.gbff.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in $(eval echo "$INPUT_DATA_PATH/*.gbk.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in $(eval echo "$INPUT_DATA_PATH/*.embl.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in $(eval echo "$INPUT_DATA_PATH/*.fasta.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in $(eval echo "$INPUT_DATA_PATH/*.fna.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in $(eval echo "$INPUT_DATA_PATH/*.fa.gz"); do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done


# create input file list

for f in $(eval echo "$INPUT_DATA_PATH/*.gbff"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done

for f in $(eval echo "$INPUT_DATA_PATH/*.gbk"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done

for f in $(eval echo "$INPUT_DATA_PATH/*.embl"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done

for f in $(eval echo "$INPUT_DATA_PATH/*.fasta"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done 

for f in $(eval echo "$INPUT_DATA_PATH/*.fna"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done

for f in $(eval echo "$INPUT_DATA_PATH/*.fa"); do
	[ -e "$f" ] && ls $f >> genomes.list && echo "found $f, adding to input list" 
done

# run mult_antiSMASH.pl
perl mult_antiSMASH.pl -i genomes.list -c $NPROC

# moves output antiSMASH directories to ../Data/antiSMASH_annotations/
if [ -d $OUTPUT_DIR ]
then
	echo "$OUTPUT_DIR found - moving antiSMASH output directories there"
else
	echo "Making $OUTPUT_DIR - moving antiSMASH directories there"
	mkdir $OUTPUT_DIR
fi
mv */ $OUTPUT_DIR
mv files.list $OUTPUT_DIR

# runs cluster_renamer.pl on newly annotated gbks
ls ${OUTPUT_DIR}/*/*cluster*.gbk > BGCs.list
perl cluster_renamer.pl -i BGCs.list -j ${OUTPUT_DIR}/files.list -c $NPROC

# moves cluster_renamer.pl output files
echo "Moving renamed clusters to $OUTPUT_DIR/BGC_gbks/"
mv BGC_gbks $OUTPUT_DIR
if [ -d "${OUTPUT_RESULTS_DIR}" ]
then
	echo "${OUTPUT_RESULTS_DIR}/ found - moving nodes.tsv and BGC_file_lookup.tsv there"
else
	echo "Making ${OUTPUT_RESULTS_DIR} - moving nodes.tsv and BGC_file_lookup.tsv there"
	mkdir ${OUTPUT_RESULTS_DIR}
fi
mv nodes.tsv ${OUTPUT_RESULTS_DIR}
mv BGC_file_lookup.tsv ${OUTPUT_RESULTS_DIR}

# runs clusters_per_genome.pl on new nodes.tsv file
perl clusters_per_genome.pl -i ${OUTPUT_RESULTS_DIR}/nodes.tsv

# moves clusters_per_genome.pl output files
echo "Moving clusters_per_genome.tsv to ${OUTPUT_RESULTS_DIR}"
mv clusters_per_genome.tsv ${OUTPUT_RESULTS_DIR}

# Cleans up src/
rm genomes.list
rm BGCs.list

# Add note to CHANGELOG.txt
echo -e "* antiSMASH_annotations.sh v3.3 finished `date +%Y-%m-%d:%H:%M:%S`\n" >> $LOG_FILE

