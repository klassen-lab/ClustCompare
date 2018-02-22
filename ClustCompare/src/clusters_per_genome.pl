#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0 - January 4, 2018 - complete rewrite, derives new table from cluster_renamer.pl output table
# v2.1 - February 21, 2018 - slight tweak to accommodate modifed nodes.tsv
#
# script counts the number of each cluster type within each genome (i.e., unique DESCRIPTION field)

use strict;
use warnings;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "clusters_per_genome.pl

DEPENDANCIES: none

USAGE:
-h	displays this usage statement (also using --help)
-i	input node table produced by cluster_renamer.pl			REQUIRED				e.g., perl cluster_renamer.pl -i nodes.tsv
-o	output tab-delimited table of cluster types in each genome	DEFAULT: clusters_per_genome.tsv	e.g., perl cluster_renamer.pl -i nodes.tsv -o genome_clusters.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N				e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder using the same name that contains the antiSMASH annotation for that genome
";

# input arguements

my %options = ();
GetOptions(
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

unless ($options{outfile}){ $options{outfile} = "clusters_per_genome.tsv"};
unless ($options{quiet}){   $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
clusters_per_genome.pl	Jonathan Klassen	v2.1	Jan 4, 2018

parameters used:
	input file = $options{infile}
	output file = $options{outfile}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Counts how many of each cluster type per genome (i.e., unique description)
#############################################################################

my %clusters; # counts each type in each genome
my %types;    # lookup of all types found
open (INFILE, $options{infile}) or die $usage;
<INFILE>; # skip header line
while (<INFILE>){
	s/\s+$//;
	my @line = split /\t/;
	my $type = $line[2];
	my $descr = $line[1];
	$clusters{$descr}{$type}++;
	$types{$type} = "y";
}

#############################################################################
# Makes output table
#############################################################################

open (OUTFILE, ">$options{outfile}") or die "Cannot open output table $options{outfile}";

# print header

print OUTFILE "Genome_description";
foreach (sort keys %types){ print OUTFILE "\t", $_ }
print OUTFILE "\n";

# print output table

foreach my $descr (sort keys %clusters){
	print OUTFILE $descr;
	foreach my $type (sort keys %types){		
		if ($clusters{$descr}{$type}){ print OUTFILE "\t$clusters{$descr}{$type}" }
		else { print OUTFILE "\t0" }
	}
	print OUTFILE "\n";
}
