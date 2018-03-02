#!/usr/bin/bash

# script to compare gene clusters by their orthologous Pfam domains
# Assumes that genomes are in ../Data/genomes and clusters are in ../Data/antiSMASH_annotations/
# Will put table in ../Data/antiSMASH_annotations and modify ../CHANGELOG.txt
# v1.0 June 22, 2017
# v2.0 January 6, 2018 - Rewrote to handle new script formats
# v3.0 February 27, 2018 - Handles parameter passing, no outputs in cwd

usage()
{
	echo "USAGE:"
	echo "Syntax: ./cluster_pfam_BBH_comparison.sh [OPTIONS]"
	echo "Options:"
	echo "    -c    number of CPU threads"
	echo "    -i    path to directory containing input data files"
	echo "    -o    path to directory where output data files should be placed"
	echo "    -r    path to directory where output results files should be placed"
	echo "    -l    path to log file to update"
	echo "    -p    path to directory that contains pfamscan.pl and the pfam HMMs"
	echo "    -h    print this help message"
}

# default parameters
NPROC=$(nproc)
INPUT_DATA_PATH=../Data/antiSMASH_annotations/BGC_gbks
OUTPUT_DIR=../Data/cluster_pfam_BBH_comparison
OUTPUT_RESULTS_DIR=../Results
PFAM_DIR=~/Tools/PfamScan
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
		-p )	shift
			PFAM_DIR=$1
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
echo "Looking for pfmascan and pfam HMMs here: $PFAM_DIR"
echo "Ouput data files will be here: $OUTPUT_DIR"
echo "Ouput data results files will be here: $OUTPUT_RESULTS_DIR"
echo "Updating log: $LOG_FILE"

# Add note to CHANGELOG.txt
echo -e "## `date +%Y-%m-%d:%H:%M:%S` \n\n* cluster_pfam_BBH_comparison.sh v3.0 has started" >> $LOG_FILE 

# check for $INPUT_DATA_PATH directory
if [ -d $INPUT_DATA_PATH ]; then
	echo "$INPUT_DATA_PATH found"
else
	echo "$INPUT_DATA_PATH not found - exiting antiSMASH_annotation.sh"
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

# check for $PFAM_DIR directory
if [ -d "$HOME/Tools/PfamScan/" ]; then
	echo "$PFAM_DIR found"
else
	echo "$PFAM_DIR not found - exiting"
	exit
fi

# create input file list of BGC gbks files
for f in $INPUT_DATA_PATH/*.gbk; do
	[ -e "$f" ] && ls $f >> $OUTPUT_DIR/cluster_gbks.list && echo "found $f, adding to input list"
done

# run mult_gbk_to_faa.pl
mkdir $OUTPUT_DIR/BGC_faas
perl mult_gbk_to_faa.pl -i $OUTPUT_DIR/cluster_gbks.list -c $NPROC -d $OUTPUT_DIR/BGC_faas

# create input file list of BGC faa files
for f in $OUTPUT_DIR/BGC_faas/*.faa; do
	[ -e "$f" ] && ls $f >> $OUTPUT_DIR/cluster_faas.list
done

# run mult_pfamscan.pl
mkdir $OUTPUT_DIR/PfamScan_output
perl mult_pfamscan.pl -i $OUTPUT_DIR/cluster_faas.list -c $NPROC -p $PFAM_DIR -d $OUTPUT_DIR/PfamScan_output

# create input file list of BGC pfamscan outputs
for f in $OUTPUT_DIR/PfamScan_output/*.pfamscan; do
	[ -e "$f" ] && ls $f >> $OUTPUT_DIR/cluster_pfamscans.list
done

# run mult_pfamscan_parser.pl
mkdir $OUTPUT_DIR/BGC_domain_faas
perl mult_pfamscan_parser.pl -i $OUTPUT_DIR/cluster_pfamscans.list -j $OUTPUT_DIR/cluster_faas.list -c 4 -d $OUTPUT_DIR/BGC_domain_faas -o $OUTPUT_RESULTS_DIR/cluster_domains.tsv

# make BLAST database of all domain faas
ls $OUTPUT_DIR/BGC_domain_faas/*.faa > $OUTPUT_DIR/cluster_domain_faas.list
cat $OUTPUT_DIR/BGC_domain_faas/*.faa > $OUTPUT_DIR/all_domains.faa
makeblastdb -dbtype prot -in $OUTPUT_DIR/all_domains.faa

# runs mult_blastp.pl
mkdir $OUTPUT_DIR/BGC_BLASTps
perl mult_blastp.pl -b $OUTPUT_DIR/all_domains.faa -i $OUTPUT_DIR/cluster_domain_faas.list -c $NPROC -d $OUTPUT_DIR/BGC_BLASTps

# create input file list of BGC faa files
for f in $OUTPUT_DIR/BGC_BLASTps/*.blastp; do
	[ -e "$f" ] && ls $f >> $OUTPUT_DIR/cluster_blastps.list
done

# runs blastp_to_BBH_list.pl
perl blastp_to_BBH_list.pl -i $OUTPUT_DIR/cluster_blastps.list -a 70 -b 70 -o $OUTPUT_RESULTS_DIR/domain_BBHs.tsv

# run and parse BBH tables
perl make_cluster_similarities.pl -i $OUTPUT_RESULTS_DIR/domain_BBHs.tsv -j $OUTPUT_RESULTS_DIR/cluster_domains.tsv -o $OUTPUT_RESULTS_DIR/raw_edges.tsv
perl cluster_similarity_table_parser.pl -a 0.3 -b 70 -c 2 -d 50 -i $OUTPUT_RESULTS_DIR/raw_edges.tsv -o $OUTPUT_RESULTS_DIR/edges.tsv

# Cleans up 
rm $OUTPUT_DIR/cluster_gbks.list
rm $OUTPUT_DIR/cluster_faas.list
rm $OUTPUT_DIR/cluster_pfamscans.list
rm $OUTPUT_DIR/all_domains.faa*
rm $OUTPUT_DIR/cluster_domain_faas.list
rm $OUTPUT_DIR/cluster_blastps.list

# Add note to CHANGELOG.txt
echo -e "* cluster_pfam_BBH_comparison.sh v3.0 finished `date +%Y-%m-%d:%H:%M:%S`\n" >> $LOG_FILE
 

