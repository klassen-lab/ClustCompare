#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - March 3/17
# v1.1 - June 22/17 - captures hybrid cluster types, removes extra row from output table
# v1.2 - June 30/17 - assumes new table each time
# v2.0 - January 2, 2018 - formal input parameters, more generalizable input
# v3.0 - February 21, 2018 - separates node table metadata from file lookups as a separate table
# v3.1 - February 28, 2018 - make temp files in outdir specified by -d instead of cwd
# v3.2 - March 2, 2018 - Fix off by one error in file names
# v3.3 - March 8, 2018 - Add knownclusterblast results
# v3.4 - March 12, 2018 - Fix thread collisions when writing output files
#
# script collates and renames the clusters from each antiSMASH genomer
# input is a list of antiSMASH clusters
# output is table listing the sources, definitions, and types of each cluster, and their new names
# also produces a single folder holding all clusters

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;
use Bio::SeqIO;

############################################################################
# Processes input arguments
############################################################################

my $usage = "cluster_renamer.pl
DEPENDANCIES: Perl \"Parallel::ForkManager\" module, BioPerl
USAGE:
-c	number of cores to use								DEFAULT: 1			e.g., perl cluster_renamer.pl -i inlist -c 4
-d	output directory for renamed cluster gbk files					DEFAULT: BGC_gbks		e.g., perl cluster_renamer.pl -i inlist -d renamed_clusters
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed			REQUIRED			e.g., perl cluster_renamer.pl -i inlist
-j	input table of source data files and corresponding antiSMASH file directories	REQUIRED			e.g., perl cluster_renamer.pl -i inlist -j ../Data/antiSMASH_annotations/files.list
-l	output tab-delimited look-up table linking clusters to source data		DEFAULT: BGC_file_lookup.tsv	e.g., perl cluster_renamer.pl -i inlist -l lookup
-o	output tab-delimited table of cluster names (node table)			DEFAULT: nodes.tsv		e.g., perl cluster_renamer.pl -i inlist -o clusters.tsv
-q	run quietly, i.e., no STDOUT (Y or N)						DEFAULT: N			e.g., perl cluster_renamer.pl -i inlist -q Y
OUTPUT FILES:
	A table specified by -o lists metadata for each cluster, such as can be used to annotate nodes in a cluster similarity network. A second table specified by -l produces a look-up table linking each cluster to their cognate genomes and scaffolds.
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"d=s"  => \$options{outdir},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"j=s"  => \$options{indirtable},
	"l=s"  => \$options{outlookup},
	"o=s"  => \$options{outnodes},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input list of paths to the cluster gbk files to be analyzed is not specified:\n, $usage" if (!$options{infile});
die "Input table of source data files and corresponding antiSMASH file directories is not specified:\n, $usage" if (!$options{indirtable});

# defaut arguments

unless ($options{cores}){      $options{cores} = 1};
unless ($options{outdir}){     $options{outdir} = "BGC_gbks"};
unless ($options{outnodes}){   $options{outnodes} = "nodes.tsv"};
unless ($options{outlookup}){  $options{outlookup} = "BGC_file_lookup.tsv"};
unless ($options{quiet}){      $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
cluster_renamer.pl	Jonathan Klassen	v3.4	Mar 12, 2018
parameters used:
	input file = $options{infile}
	input table listing annotation directories = $options{indirtable}
	output directory for renamed BGC gbks = $options{outdir}\/
	output nodes table = $options{outnodes}
	output file lookup table = $options{outlookup}
	number of cores = $options{cores}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @clusters;
open (INLIST, $options{infile}) or die "Cannot open $options{infile}\n$usage";

while (<INLIST>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @clusters, $_;
}

#############################################################################
# Loads table linking genomes to antiSMASH output directories
#############################################################################

my %clusterdirs;
unless ($options{indirtable} eq "Not specified"){
	open (INDIRS, "$options{indirtable}") or die "Cannot open $options{indirtable}\n$usage";
	<INDIRS>; # assumes header line
	while (<INDIRS>){

		s/\s+$//;
		my @line = split /\t/, $_;
		$clusterdirs{$line[1]} = $line[0];
	}
}

#############################################################################
# Multithreaded cluster parsing and output file creation
#############################################################################

# check if output directory exists

if (-d $options{outdir}){ die "$options{outdir}\/ already exists, existing cluster_renamer.pl\n";}

# else make new output directory

mkdir $options{outdir} or die "Cannot make $options{outdir}";

# multithreaded cluster parsing

my $pm = new Parallel::ForkManager($options{cores});
sub cluster_parser($);
for my $a (0..$#clusters){
	$pm->start and next;
	my @return = cluster_parser($a);
	my $num = $a + 1;
	system "cp $return[2] $options{outdir}/cluster$num.gbk";
	$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# Collate data from temp output file
#############################################################################

my %collated_data;
my @node_table_header = ("Source_description", "Cluster_file", "Cluster_type", "Source_contig_accession"); 
for my $a (0..$#clusters){
	my $num = $a + 1;

	open (INTEMP, "$options{outdir}/$num.temp") or die "Cannot open $options{outdir}/$num.temp";

	while (<INTEMP>){
		s/\s+$//;
		s/1\. (BGC\w*)/$1/;
		s/_biosyn\w*/_BGC\t/;          # adjust knowncluster column
		s/\((.*)\%.*$/$1/;
	#	s/\%.*/\%/;
		my @line = split /\t/, $_;
		$collated_data{$line[0]}{Source_description} = $line[1];
		$collated_data{$line[0]}{Cluster_file} = $line[2];
		$collated_data{$line[0]}{Cluster_type} = $line[3];
		$collated_data{$line[0]}{Source_contig_accession} = $line[4];
		if ($line[5]){
			$node_table_header[4] = "Top_known_cluster_id";
			$collated_data{$line[0]}{Top_known_cluster_id} = $line[5];
			$node_table_header[5] = "Top_known_cluster_name";
			$collated_data{$line[0]}{Top_known_cluster_name} = $line[6];
			$node_table_header[6] = "Top_known_cluster_score";
			$collated_data{$line[0]}{Top_known_cluster_score} = $line[7];
		}
	}
	close INTEMP;
	system "rm $options{outdir}/$num.temp";
}

#############################################################################
# Generate node output file
#############################################################################

open (OUTNODES, ">$options{outnodes}") or die "Cannot open output node file $options{outnodes}\n$usage";
print OUTNODES "Cluster_id\t";
foreach my $node_table_key (@node_table_header){
	print OUTNODES "$node_table_key\t";
}
print OUTNODES "\n";
#print OUTNODES "Cluster_id\tSource_description\tCluster_type\tTop_known_cluster_id\tTop_known_cluster_name\tTop_known_cluster_score\n";
foreach my $cluster_id (sort {$a <=> $b} keys %collated_data){
#	my $num = $cluster_id + 1;
	print OUTNODES "cluster$cluster_id";
	foreach my $node_table_key (@node_table_header){
		print OUTNODES "\t$collated_data{$cluster_id}{$node_table_key}";
	}
	print OUTNODES "\n";
}

# generate lookup table output file

open (OUTLOOKUP, ">$options{outlookup}") or die "Cannot open output lookup table file $options{outlookup}\n$usage";
print OUTLOOKUP "Cluster_id\tCluster_file\tSource_genome_file\tSource_contig_accession\n";
foreach my $cluster_id (sort {$a <=> $b} keys %collated_data){
	#my $num = $cluster_id + 1;
	print OUTLOOKUP "cluster$cluster_id\t$collated_data{$cluster_id}{Cluster_file}\t";
	my $dirname = $collated_data{$cluster_id}{Cluster_file};
	$dirname =~ s/\/(\w|\.|\[|\]|\(|\))+$//;
	$dirname =~ s/^.+\///;
	if ($clusterdirs{$dirname}){
		print OUTLOOKUP "$clusterdirs{$dirname}\t";
	}
	else {
		print OUTLOOKUP "NA\t";
	}
	print OUTLOOKUP "$collated_data{$cluster_id}{Source_contig_accession}\n";
}

#############################################################################
# Subroutine to collate BGCs 
#############################################################################

sub cluster_parser($){
	my $counter = $_[0];								# cluster id
	my $num = $counter + 1;

	my @return;
	if ($options{quiet} eq "N"){ print "Collating cluster $clusters[$counter]: ", $counter + 1, " of ", scalar @clusters, "\n"};
	my $filepath = $clusters[$counter]; 						# filepath

	my ($descr, $type, $accession, $top_known);
	my $seqio_obj = Bio::SeqIO->new(-file => $clusters[$counter], -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # loop through contigs
		$descr = $seq_obj->desc; 						# contig description
		$accession = $seq_obj->accession_number;			 	# contig accession number
		for my $feat_obj ($seq_obj->get_SeqFeatures){ # loop through annotation features
			if ($feat_obj->primary_tag eq "cluster"){       
				my @cluster_type; 
				my @knownclusters;
				push @cluster_type, $feat_obj->get_tag_values("product");
				push @knownclusters, $feat_obj->get_tag_values("knownclusterblast")
				    if ($feat_obj->has_tag("knownclusterblast"));
				$top_known = $knownclusters[0];
				$type = $cluster_type[0];				# BGC type
			}
		}	
	}

	# collates summary lines as temp files to avoid thread collisions

	open (OUTSUMMARY, ">$options{outdir}/$num.temp") or die "Cannot open $options{outdir}/$num.temp";
	print OUTSUMMARY "$num\t$descr\t$filepath\t$type\t$accession\t";
		if ($top_known){
			print OUTSUMMARY "$top_known\n";
		}
		else {
			print OUTSUMMARY "\n";
		}
	close OUTSUMMARY;
	
	return ($counter, $descr, $filepath, $type, $accession, $top_known); 
}
