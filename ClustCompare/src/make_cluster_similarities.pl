#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - June 29/17
# v2.0 - January 8, 2017 - general input/outputs, formal parameter passing. No longer adds cluster types
#
# Parses output for blastp_to_BBH.pl to determine similarities between each cluster. 
# Gets domain info from output of mult_pfamscan_parser.pl


use strict;
use warnings;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "make_cluster_similarities.pl

DEPENDANCIES: none

USAGE:
-h	displays this usage statement (also using --help)
-i	input BBH table output by blastp_to_BBH.pl						REQUIRED		e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv
-j	input table of domains per cluster output by mult_pfamscan_parser.pl			REQUIRED		e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv
-o	output tab-delimited table of domain orthologs shared by two clusters (network edges)	DEFAULT: raw_edges.tsv	e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv -o BBHs.tsv
-q	run quietly, i.e., no STDOUT (Y or N)							DEFAULT: N		e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv -q Y

OUTPUT FILES:
	a table listing the number of domains that are shared between two BGCs (out of all possible domains) and how closely these domains are related to each other; this table contains all possible relationships without filtering
";

# input arguements

my %options = ();
GetOptions( 
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{inBBH},
	"j=s"  => \$options{indomains},
	"o=s"  => \$options{outfile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input BBH table is not specified:\n, $usage" if (!$options{inBBH});
die "Input domains per cluster table is not specified:\n, $usage" if (!$options{indomains});

# defaut arguments

unless ($options{outfile}){ $options{outfile} = "raw_edges.tsv"};
unless ($options{quiet}){   $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
make_cluster_similarities.pl	Jonathan Klassen	v2.0	Jan 8, 2018

parameters used:
	input BBH table = $options{inBBH}
	input cluster domain table = $options{indomains}
	output file = $options{outfile}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @bbhs;
open (INBBHs, $options{inBBH}) or die $usage;
while (<INBBHs>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @bbhs, $_;
}
shift @bbhs; # skip header

my @domains;
open (INDOMAINS, $options{indomains}) or die $usage;
while (<INDOMAINS>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @domains, $_;
}
shift @domains; # skip header

#############################################################################
# Makes lookup table counting the number of domains in each cluster
#############################################################################

my %num_cluster_domains;
foreach my $cluster (@domains){
	$cluster =~ s/\s*$//;
	my @line = split /\t/, $cluster;
	$num_cluster_domains{$line[0]} = $line[1];
}

#############################################################################
# Counts number of bbh hits between each pair of clusters
#############################################################################

my %cluster_similarities;
foreach my $bbh_line (@bbhs){
	$bbh_line =~ s/\s*$//;
	my @line = split /\t/, $bbh_line;
	$line[0] =~ s/_\w+$//;
	$line[1] =~ s/_\w+$//;
	if ($cluster_similarities{$line[1]}{$line[0]}){
		push @{$cluster_similarities{$line[1]}{$line[0]}}, $line[2];
	}
	else {
		push @{$cluster_similarities{$line[0]}{$line[1]}}, $line[2];
	}
}

#############################################################################
# Create output table
#############################################################################

open (OUTFILE, ">$options{outfile}") or die "Cannot open $options{outfile}";
print OUTFILE "Cluster_1\tCluster_2\tScore\tAvg_domain_identity\tNum_shared_domains\tPercent_shared_domains\tDomains_in_cluster_1\t",
	"Domains_in_cluster_2\tPercent_ID_shared_domains\n";
foreach my $cluster1 (sort keys %cluster_similarities){
	foreach my $cluster2 (sort keys %{$cluster_similarities{$cluster1}}){
		my $avg_id;
		foreach (@{$cluster_similarities{$cluster1}{$cluster2}}){
			$avg_id += $_;
		}
		my $num_shared_domains = scalar @{$cluster_similarities{$cluster1}{$cluster2}};
		$avg_id = $avg_id / $num_shared_domains;
		my $min_domains;
		if ($num_cluster_domains{$cluster1} <= $num_cluster_domains{$cluster2}){ 
			$min_domains = $num_cluster_domains{$cluster1};
		}
		else {
			$min_domains = $num_cluster_domains{$cluster2};
		}
		my $percent_shared_domains = $num_shared_domains / $min_domains * 100;
		my $score = ($avg_id / 100) * ($percent_shared_domains / 100);
		print OUTFILE "$cluster1\t$cluster2\t$score\t$avg_id\t$num_shared_domains\t$percent_shared_domains\t",
			"$num_cluster_domains{$cluster1}\t$num_cluster_domains{$cluster2}\t",
			join ",", @{$cluster_similarities{$cluster1}{$cluster2}}, "\n";
	}
}



		



