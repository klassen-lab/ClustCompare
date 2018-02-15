#!/usr/bin/bash

# Script to annotate a batch of NCBI genomes using antiSMASH. Assumes that run_antismash.py works from terminal
# Will put annotations in ../Data/ and modify ../CHANGELOG.txt
# v3.0 January 2, 2018 - adapted to Jan 2018 script versions. Requires all genomes be located in ../Data/genomes. Runs on .gbk, .gbff, .fasta, .fna, .fa nucleotide formats.
# v3.1 February 15, 2018 - Accommodates genomes.list output file

# Add note to CHANGELOG.txt
echo -e "## `date +%Y-%m-%d:%H:%M:%S` \n\n* antiSMASH_annotations.sh v3.0 started" >> ../CHANGELOG.txt 

# check for ../Data/genomes directory
if [ -d "../Data/genomes/" ]; then
	echo "../Data/genomes found"
else
	echo "../Data/genomes not found - exiting antiSMASH_annotation.sh"
	exit
fi

# find number of cores
NPROC=$(nproc)
echo "Using $NPROC cores"

# unzip any zipped files, if they exist
for f in ../Data/genomes/*.gbff.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in ../Data/genomes/*.gbk.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in ../Data/genomes/*.embl.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in ../Data/genomes/*.fasta.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in ../Data/genomes/*.fna.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done

for f in ../Data/genomes/*.fa.gz; do
	[ -e "$f" ] && gunzip $f && echo "unzipping $f"
done


# create input file list

for f in ../Data/genomes/*.gbff; do
	[ -e "$f" ] && ls $f >> genomes.list
done

for f in ../Data/genomes/*.gbk; do
	[ -e "$f" ] && ls $f >> genomes.list
done

for f in ../Data/genomes/*.embl; do
	[ -e "$f" ] && ls $f >> genomes.list
done

for f in ../Data/genomes/*.fasta; do
	[ -e "$f" ] && ls $f >> genomes.list
done

for f in ../Data/genomes/*.fna; do
	[ -e "$f" ] && ls $f >> genomes.list
done

for f in ../Data/genomes/*.fa; do
	[ -e "$f" ] && ls $f >> genomes.list
done

# run mult_antiSMASH.pl
perl mult_antiSMASH.pl -i genomes.list -c $NPROC

# moves output antiSMASH directories to ../Data/antiSMASH_annotations/
if [ -d "../Data/antiSMASH_annotations/" ]
then
	echo "../Data/antiSMASH_annotations found - moving antiSMASH output directories there"
else
	echo "Making ../Data/antiSMASH_annotations - moving antiSMASH directories there"
	mkdir ../Data/antiSMASH_annotations
fi
mv */ ../Data/antiSMASH_annotations
mv files.list ../Data/antiSMASH_annotations

# runs cluster_renamer.pl on newly annotated gbks
ls ../Data/antiSMASH_annotations/*/*cluster*.gbk > BGCs.list
perl cluster_renamer.pl -i BGCs.list -c $NPROC

# moves cluster_renamer.pl output files
echo "Moving renamed clusters to ../Data/BGC_gbks/"
mv BGC_gbks ../Data/
if [ -d "../Results/" ]
then
	echo "../Results/ found - moving nodes.tsv there"
else
	echo "Making ../Results - moving nodes.tsv there"
	mkdir ../Results
fi
mv nodes.tsv ../Results/

# runs clusters_per_genome.pl on new nodes.tsv file
perl clusters_per_genome.pl -i ../Results/nodes.tsv

# moves clusters_per_genome.pl output files
echo "Moving clusters_per_genome.tsv to ../Results"
mv clusters_per_genome.tsv ../Results

# Cleans up src/
rm genomes.list
rm BGCs.list

# Add note to CHANGELOG.txt
echo -e "* antiSMASH_annotations.sh v3.0 finished `date +%Y-%m-%d:%H:%M:%S`\n" >> ../CHANGELOG.txt
 

