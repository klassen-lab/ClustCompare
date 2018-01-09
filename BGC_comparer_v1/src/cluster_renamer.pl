#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - March 3/17
# v1.1 - June 22/17 - captures hybrid cluster types, removes extra row from output table
# v1.2 - June 30/17 - assumes new table each time
# v2.0 - January 2, 2018 - formal input parameters, more generalizable input
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
-c	number of cores to use						DEFAULT: 1		e.g., perl cluster_renamer.pl -i inlist -c 4
-d	output directory for renamed cluster gbk files			DEFAULT: BGC_gbks	e.g., perl cluster_renamer.pl -i inlist -d renamed_clusters
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl cluster_renamer.pl -i inlist
-o	output tab-delimited table of cluster names (node table)	DEFAULT: nodes.tsv	e.g., perl cluster_renamer.pl -i inlist -o clusters.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder using the same name that contains the antiSMASH annotation for that genome
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"d=s"  => \$options{outdir},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"o=s"  => \$options{outnodes},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n, $usage" if (!$options{infile});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outdir}){   $options{outdir} = "BGC_gbks"};
unless ($options{outnodes}){ $options{outnodes} = "nodes.tsv"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
cluster_renamer.pl	Jonathan Klassen	v2.0	Jan 2, 2018

parameters used:
	input file = $options{infile}
	output directory for renamed BGC gbks = $options{outdir}\/
	number of cores = $options{cores}
	output nodes table = $options{outnodes}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @clusters;
open (INLIST, $options{infile}) or die $usage;

while (<INLIST>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @clusters, $_;
}

#############################################################################
# Multithreaded cluster parsing and output file creation
#############################################################################

# check if output directory exists

if (-d $options{outdir}){ die "$options{outdir}\/ already exists, existing cluster_renamer.pl\n";}

# else make new output directory

mkdir $options{outdir} or die "Cannot make $options{outdir}";

# open output file

open (OUTNODES, ">$options{outnodes}") or die "Cannot open output node file $options{outnodes}";
print OUTNODES "Cluster_id\tSource_description\tSource_file\tCluster_type\tCluster_contig_accession\n";

# multithreaded cluster parsing

my $pm = new Parallel::ForkManager($options{cores});
sub cluster_parser($);
for my $a (0..$#clusters){
	$pm->start and next;
	my @return = cluster_parser($a);
	system "cp $return[2] $options{outdir}/cluster$a.gbk";
	print OUTNODES "cluster$a\t$return[1]\t$return[2]\t$return[3]\t$return[4]\n";
	$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# Subroutine to collate BGCs 
#############################################################################

sub cluster_parser($){
	my $counter = $_[0];								# cluster id
	my @return;
	if ($options{quiet} eq "N"){ print "Collating cluster $clusters[$counter]: ", $counter + 1, " of ", scalar @clusters, "\n"};
	my $filepath = $clusters[$counter]; 						# filepath

	my ($descr, $type, $accession);
	my $seqio_obj = Bio::SeqIO->new(-file => $clusters[$counter], -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # loop through contigs
		$descr = $seq_obj->desc; 						# contig description
		$accession = $seq_obj->accession_number;			 	# contig accession number
		for my $feat_obj ($seq_obj->get_SeqFeatures){ # loop through annotation features
			if ($feat_obj->primary_tag eq "cluster"){       
				my @cluster_type;
				push @cluster_type, $feat_obj->get_tag_values("product");	
				$type = $cluster_type[0];				# BGC type
			}
		}	
	}	
	return ($counter, $descr, $filepath, $type, $accession); 
}


