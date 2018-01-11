#!/usr/bin/bash

# script to compare gene clusters by their orthologous Pfam domains
# Assumes that genomes are in ../Data/genomes and clusters are in ../Data/antiSMASH_annotations/
# Will put table in ../Data/antiSMASH_annotations and modify ../CHANGELOG.txt
# v1.0 June 22, 2017
# v2.0 January 6, 2018 - Rewrote to handle new script formats

# Add note to CHANGELOG.txt
echo -e "## `date +%Y-%m-%d:%H:%M:%S` \n\n* cluster_pfam_BBH_comparison.sh v2.0 has started" >> ../CHANGELOG.txt 

# check for ../Data/genomes directory
if [ -d "../Data/BGC_gbks/" ]; then
	echo "../Data/BGC_gbks/ found"
else
	echo "../Data/BGC_gbks/ not found - exiting cluster_pfam_BBH_comparison.sh"
	exit
fi

# check for ~/Tools/PfamScan (needed for mult_pfamscan.pl)
if [ -d "$HOME/Tools/PfamScan/" ]; then
	echo "~/Tools/PfamScan/ found, assuming if contains pfamscan.pl and the pfam HMMs"
else
	echo "~/Tools/PfamScan/ not found, modify this script to contain the correct directory that contains pfamscan.pl and the pfam HMMs"
	exit
fi

# find number of cores
NPROC=$(nproc)
echo "Using $NPROC cores"

# create input file list of BGC gbks files
for f in ../Data/BGC_gbks/*.gbk; do
	[ -e "$f" ] && ls $f >> cluster_gbks.list
done

# run mult_gbk_to_faa.pl
perl mult_gbk_to_faa.pl -i cluster_gbks.list -c $NPROC

# moves folder of BGC faa files to ../Data/
echo "Moving cluster faa files to ../Data/BGC_faas/"
mv BGC_faas ../Data/

# create input file list of BGC faa files
for f in ../Data/BGC_faas/*.faa; do
	[ -e "$f" ] && ls $f >> cluster_faas.list
done

# run mult_pfamscan.pl
perl mult_pfamscan.pl -i cluster_faas.list -c $NPROC -p ~/Tools/PfamScan

# moves folder of pfamscan output files to ../Data/
echo "Moving pfamscan output files to ../Data/BGC_pfamscans/"
mv BGC_pfamscans ../Data/

# create input file list of BGC pfamscan outputs
for f in ../Data/BGC_pfamscans/*.pfamscan; do
	[ -e "$f" ] && ls $f >> cluster_pfamscans.list
done

# run mult_pfamscan_parser.pl
perl mult_pfamscan_parser.pl -i cluster_pfamscans.list -j cluster_faas.list -c 4

# moves mult_pfamscan_parser.pl output files
echo "Moving BGC_domain_faas ../Data/"
mv BGC_domain_faas ../Data/

if [ -d "../Results/" ]
then
	echo "../Results/ found - moving cluster_domains.tsv there"
else
	echo "Making ../Results - moving cluster_domains.tsv there"
	mkdir ../Results
fi
mv cluster_domains.tsv ../Results/

# make BLAST database of all domain faas
cat ../Data/BGC_domain_faas/*.faa > all_domains.faa
makeblastdb -dbtype prot -in all_domains.faa
ls ../Data/BGC_domain_faas/*.faa > cluster_domain_faas.list

# runs mult_blastp.pl
perl mult_blastp.pl -b all_domains.faa -i cluster_domain_faas.list -c $NPROC

# moves folder of BLASTp output files to ../Data/
echo "Moving BGC_BLASTps/ to ../Data/"
mv BGC_BLASTps ../Data/

# create input file list of BGC faa files
for f in ../Data/BGC_BLASTps/*.blastp; do
	[ -e "$f" ] && ls $f >> cluster_blastps.list
done

# runs blastp_to_BBH_list.pl
perl blastp_to_BBH_list.pl -i cluster_blastps.list -a 70 -b 70

# moves BBH output table to ../Results
echo "Moving domain_BBHs.tsv to ../Results/"
mv domain_BBHs.tsv ../Results/

# run and parse BBH tables
perl make_cluster_similarities.pl -i ../Results/domain_BBHs.tsv  -j ../Results/cluster_domains.tsv 
echo "Moving raw_edges.tsv to ../Results/"
mv raw_edges.tsv ../Results/
perl cluster_similarity_table_parser.pl -a 0.3 -b 70 -c 2 -d 50 -i ../Results/raw_edges.tsv
echo "Moving edges.tsv to ../Results/"
mv edges.tsv ../Results/

# Cleans up src/
rm cluster_gbks.list
rm cluster_faas.list
rm cluster_pfamscans.list
rm all_domains.faa*
rm cluster_domain_faas.list
rm cluster_blastps.list

# Add note to CHANGELOG.txt
echo -e "* cluster_pfam_BBH_comparison.sh v2.0 finished `date +%Y-%m-%d:%H:%M:%S`\n" >> ../CHANGELOG.txt
 

