#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - June 29/17
# v2.0 - January 8/17 - general input/output, formal input parameters
#
# Parses ../Data/antiSMASH_annotations/domain_comparisons/cluster_similarities.table to exclude low-quality matches

use strict;
use warnings;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "clusters_similarly_table_parser.pl

DEPENDANCIES: none

USAGE:
-a	minimum cluster similarity score threshold (range: 0-1)				DEFAULT: 0.3
-b	minimum average ortholog percent id threshold					DEFAULT: 70
-c	minimum number of shared domains between clusters				DEFAULT: 2
-d	minimum percent of shared domains between clusters				DEFAULT: 50

-h	displays this usage statement (also using --help)
-i	input cluster similariy table produced by make_cluster_similarities.pl		REQUIRED		e.g., perl clusters_similarly_table_parser.pl -i raw_edges.tsv
-o	output tab-delimited table of parsed cluster similarities (network edges)	DEFAULT: edges.tsv	e.g., perl clusters_similarly_table_parser.pl -i raw_edges.tsv -o cluster_edges.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N				e.g., perl clusters_similarly_table_parser.pl -i inlist -q Y

OUTPUT FILES:
	a table listing the number of domains that are shared between two BGCs (out of all possible domains) and how closely these domains are related to each other; this table contains only relationships passing the specified threshold filters
";

# input arguements

my %options = ();
GetOptions(
	"a=s"  => \$options{score},
        "b=i"  => \$options{id},
        "c=i"  => \$options{shared},
        "d=i"  => \$options{percent_doms},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"o=s"  => \$options{outfile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n, $usage" if (!$options{infile});

# defaut arguments

unless ($options{score}){        $options{score} = 0.3};
unless ($options{id}){           $options{id} = 70};
unless ($options{shared}){       $options{shared} = 2};
unless ($options{percent_doms}){ $options{percent_doms} = 50};
unless ($options{outfile}){      $options{outfile} = "edges.tsv"};
unless ($options{quiet}){        $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
clusters_similarly_table_parser.pl	Jonathan Klassen	v2.0	Jan 4, 2018

parameters used:
	minimum score threshold = $options{score}
	minimum average ortholog % identity = $options{id}
	minimum number of shared domains = $options{shared}
	minimum % domains shared = $options{percent_doms}
	input file = $options{infile}
	output file = $options{outfile}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");


#############################################################################
# Loads list of input files
#############################################################################

my @infile;
open (INFILE, $options{infile}) or die $usage;
while (<INFILE>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @infile, $_;
}

#############################################################################
# Create output file
#############################################################################

open (OUTFILE, ">$options{outfile}") or die "Cannot open $options{outfile}";
print OUTFILE "$infile[0]\n";
shift @infile; # skip header

# only include selected lines of input file meeting threshold cutoffs
foreach (@infile){
	my @line = split /\t/;
	next if ($line[2] < $options{score});
	next if ($line[3] < $options{id});
	next if ($line[4] < $options{shared});
	next if ($line[5] < $options{percent_doms});
	print OUTFILE $_, "\n";
}

