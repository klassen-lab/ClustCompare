#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - June 29/17
# v2.0 - January 8, 2017 - improved input parameter parsing, more generalized input/output
#
# Parses multiple BLASTp outputs against a common database to generate a list of bidirectional best BLAST hits
# Stores output table in ../Data/antiSMASH_annotations/domain_comparisons/BLASTp_vs_all_clusters/

use strict;
use warnings;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "blastp_to_BBH_list.pl

DEPENDANCIES: none

USAGE:
-a	minimum percent identify to declare a BBH match						DEFAULT: 70			e.g., perl blastp_to_BBH_list.pl -i inlist -a 50
-b	minimum percent sequence overlap to declare a BBH match					DEFAULT: 70			e.g., perl blastp_to_BBH_list.pl -i inlist -b 80
-h	displays this usage statement (also using --help)
-i	input list of blastp output files (modified outfmt 7) produced by mult_blastp.pl	REQUIRED			e.g., perl blastp_to_BBH_list.pl -i inlist
-o	output tab-delimited table of cluster types in each genome				DEFAULT: domain_BBHs.tsv	e.g., perl blastp_to_BBH_list.pl -i inlist -o BBHs.tsv
-q	run quietly, i.e., no STDOUT (Y or N)							DEFAULT: N			e.g., perl blastp_to_BBH_list.pl -i inlist -q Y

OUTPUT FILES:
	a table listing the single best BBHs for each domain comparison between two BGCs, and their percent identify to each other
";

# input arguements

my %options = ();
GetOptions( 
	"a=i"  => \$options{id},
	"b=i"  => \$options{overlap},
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

unless ($options{id}){      $options{id} = 70};
unless ($options{overlap}){ $options{overlap} = 70};
unless ($options{outfile}){ $options{outfile} = "domain_BBHs.tsv"};
unless ($options{quiet}){   $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
blastp_to_BBH.pl	Jonathan Klassen	v2.0	Jan 8, 2018

parameters used:
	percent sequence id threshold = $options{id}
	percent sequence overlap threshold = $options{overlap}
	input file = $options{infile}
	output file = $options{outfile}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @infiles;
open (INLIST, $options{infile}) or die $usage;

while (<INLIST>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @infiles, $_;
}

############################################################################
# Extract best matches from all blastp files, stores as giant hash
############################################################################

my $count = 0;
my %best_hits;
foreach my $infile (@infiles){
	$infile =~ s/\s*$//;
	open (INFILE, $infile) or die "Cannot open $infile";
	$count++;
	print "Parsing $infile: $count of ", scalar @infiles, "\n";
	while (<INFILE>){
		next if (/^#/);
		my @line = split /\t/;
		next if (scalar @line != 14);
		next if ($line[0] eq $line[1]); # no self hits
		next if ($best_hits{$line[0]}{$line[1]} and $best_hits{$line[0]}{$line[1]} < $line[2]); # only top hits kept
		next if ($line[2] < $options{id}); # no low %ID hits
		(my $name1 = $line[0]) =~ s/_\w*$//; # no hits to own cluster
		(my $name2 = $line[1]) =~ s/_\w*$//;
		next if ($name1 eq $name2);
		my $max_length = 0; # no low % overlap hits
		if ($line[12] >= $line[13]){
			$max_length = $line[12];
		}
		else {
			$max_length = $line[13];
		}
		next if ($line[3] / $max_length * 100 < $options{overlap});
		$best_hits{$line[0]}{$line[1]} = $line[2]; # keep %ID if all thresholds passed
	}
	close INFILE;
}

############################################################################
# Creates output table of all BBHs
############################################################################

print "\nFinding BBHs\n";
open (OUTFILE, ">$options{outfile}") or die "Cannot open $options{outfile}";

print OUTFILE "Domain_1\tDomain_2\tPercent_ID\n";
my %found_bbh;
foreach my $hit1 (sort keys %best_hits){
	foreach my $hit2 (sort keys %{$best_hits{$hit1}}){
		next unless ($best_hits{$hit1}{$hit2} and $best_hits{$hit2}{$hit1}); # require bidirectional hits
		next if ($found_bbh{$hit1}{$hit2} or $found_bbh{$hit2}{$hit1});
		my $avg_id = ($best_hits{$hit1}{$hit2} + $best_hits{$hit2}{$hit1}) / 2;
		print OUTFILE "$hit1\t$hit2\t$avg_id\n";
		$found_bbh{$hit1}{$hit2} = $found_bbh{$hit2}{$hit1} = "y";
	}
}
