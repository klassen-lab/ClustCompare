#!/usr/bin/bash

# Script to determine the fragmentation of BGCs by aligning them against their source genomes
# Assumes antiSMASH_annotation.sh has been run first (or equivilent)

# v1.0 February 21, 2018 
# v1.1 February 27, 2018 - More general file handling

usage()
{
	echo "USAGE:"
	echo "Syntax: ./find_cluster_completeness.sh [OPTIONS]"
	echo "Options:"
	echo "    -c    number of CPU threads"
	echo "    -i    path to input lookup table linking BGCs to their source genomes"
	echo "    -j    path to input nodes.tsv file for updating"
	echo "    -o    path to directory where output data files should be placed"
	echo "    -r    path to directory where output results files should be placed"
	echo "    -l    path to log file to update"
	echo "    -h    print this help message"
}

# default parameters
NPROC=$(nproc)
INPUT_LOOKUP_PATH=../Results/BGC_file_lookup.tsv
INPUT_NODES_PATH=../Results/nodes.tsv
OUTPUT_DIR=../Data/BGC_fragmentation
OUTPUT_RESULTS_DIR=../Results
LOG_FILE=../CHANGELOG.txt

# parse command line options
while [ "$1" != "" ]; do
	case $1 in
		-c )	shift
			NPROC=$1
			;;
		-i )	shift
			INPUT_LOOKUP_PATH=$1
			;;
		-j )	shift
			INPUT_NODES_PATH=$1
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
echo "Looking for input lookup table here: $INPUT_LOOKUP_PATH"
echo "Looking for input nodes table here: $INPUT_NODES_PATH"
echo "Ouput data files will be here: $OUTPUT_DIR"
echo "Ouput data results files will be here: $OUTPUT_RESULTS_DIR"
echo "Updating log: $LOG_FILE"

# Add note to log
echo -e "## `date +%Y-%m-%d:%H:%M:%S` \n\n* find_cluster_completeness.sh v1.1 started" >> $LOG_FILE 

# check for $INPUT_LOOKUP_PATH file
if [ -e $INPUT_LOOKUP_PATH ]; then
	echo "$INPUT_LOOKUP_PATH found"
else
	echo "$INPUT_LOOKUP_PATH not found - find_cluster_completeness.sh"
	exit
fi

# check for $INPUT_NODES_PATH file
if [ -e $INPUT_NODES_PATH ]; then
	echo "$INPUT_NODES_PATH found"
else
	echo "$INPUT_NODES_PATH not found - find_cluster_completeness.sh"
	exit
fi

# check for $OUTPUT_DIR directory
if [ -d $OUTPUT_DIR ]
then
	echo "$OUTPUT_DIR found"
else
	echo "$OUTPUT_DIR not found - creating"
	mkdir $OUTPUT_DIR
fi

# check for $OUTPUT_RESULTS_DIR directory
if [ -d $OUTPUT_RESULTS_DIR ]
then
	echo "$OUTPUT_RESULTS_DIR found"
else
	echo "$OUTPUT_RESULTS_DIR not found - creating"
	mkdir $OUTPUT_RESULTS_DIR
fi
# run find_cluster_completeness.pl
perl find_cluster_completeness.pl -i $INPUT_LOOKUP_PATH -j $INPUT_NODES_PATH -c $NPROC -n $OUTPUT_RESULTS_DIR/nodes.tsv -o $OUTPUT_RESULTS_DIR/BGC_fragmentation.tsv -p $OUTPUT_DIR

# Add note to CHANGELOG.txt
echo -e "* find_cluster_completeness.sh v1.1 finished `date +%Y-%m-%d:%H:%M:%S`\n" >> $LOG_FILE

