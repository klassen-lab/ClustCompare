#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0 - January 6, 2018 - formal input parameters, more generalizable input, multithreaded parsing of pfamscan outputs (not gathering faas)
# v2.1 - February 27, 2018 - Slight modification to outdir creation
#
# script uses parses a list of PfamScan output files to create a summary table of domains in each input BGC
# also extracts the faa sequence for each detected domain for comparison to each other in subsequent steps

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "mult_pfamscan_parser.pl

DEPENDANCIES: Perl \"Parallel::ForkManager\" module

USAGE:
-c	number of cores to use						DEFAULT: 1			e.g., perl mult_pfamscan.pl -i inlist -j infaas c 4
-d	output directory for the pfam domain faas for each cluster	DEFAULT: BGC_domain_faas	e.g., perl mult_pfamscan.pl -i inlist -j infaas -d BGC_domains
-h	displays this usage statement (also using --help)
-i	input list of the cluster pfamscan files to be analyzed		REQUIRED			e.g., perl mult_pfamscan.pl -i inlist -j infaas
-j	input list of cluster faa files on which pfamscan was run	REQUIRED			e.g., perl mult_pfamscan.pl -i inlist -j infaas
-o	output summary file of pfam domains for each cluster		DEFAULT: cluster_domains.tsv	e.g., perl mult_pfamscan.pl -i inlist -j infaas -o domains.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N			e.g., perl mult_pfamscan.pl -i inlist -j infaas q Y

OUTPUT FILES:
	a folder containing faa files of the domains annotated by pfam in each cluster, and a table listing the domains in each cluster.
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"d=s"  => \$options{outdir},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"j=s"  => \$options{infaas},
	"o=s"  => \$options{outfile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input list of pfamscan files is not specified:\n $usage" if (!$options{infile});
die "Input list of cluster faa files is not specified:\n $usage" if (!$options{infaas});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outdir}){   $options{outdir} = "BGC_domain_faas"};
unless ($options{outfile}){  $options{outfile} = "cluster_domains.tsv"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_pfamscan_parser.pl		Jonathan Klassen	v2.1	Feb 27, 2018

parameters used:
	input list of pfamscan files = $options{infile}
	input list of faa files = $options{infaas}
	output directory for domain faa output files = $options{outdir}\/
	output table of domains in each cluster = $options{outfile}
	number of cores = $options{cores}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads lists of input files
#############################################################################

my @inpfamscans;
open (INPFAMSCAN, $options{infile}) or die $usage;

while (<INPFAMSCAN>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @inpfamscans, $_;
}

my @infaas;
open (INFAAS, $options{infaas}) or die $usage;

while (<INFAAS>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @infaas, $_;
}

#############################################################################
# Set up hash from which to extract domain faa sequences
#############################################################################

my %faas;
for my $counter (0..$#infaas){
	$infaas[$counter] =~ s/\s*$//;
	print "Collecting faa sequences from $infaas[$counter]: ", $counter + 1, " of ", scalar @infaas, "\n";
	(my $name = $infaas[$counter]) =~ s/\.\w*$//;
	$name =~ s/^.*\///; 
	open (INFAA, $infaas[$counter]) or die "Cannot open $infaas[$counter]";
	my $header;
	while (<INFAA>){
		s/\s*$//;
		if (/^>/){
			s/^>//;
			$header = $_;
		}
		else {
			$faas{$name}{$header} = $_;
		}
	}
	close INFAA;
}	

#############################################################################
# Multithreaded pfamscan
#############################################################################

# check if output directory exists

unless (-d $options{outdir}){ 
	mkdir $options{outdir} or die "Cannot make $options{outdir}";
}

# multithreaded pfamscan

my $pm = new Parallel::ForkManager($options{cores});
sub run_pfamscan_parser($);
for my $a (0..$#inpfamscans){
	$pm->start and next;
	my @return = run_pfamscan_parser($a);
		$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# Start output table
#############################################################################

open (OUTTABLE, ">$options{outfile}") or die "Cannot open $options{outfile}";
print OUTTABLE "Cluster_id\tNumber_of_domains\tDomain_ids\n";
for $counter (0..$#inpfamscans){
	open (INTEMP, "$options{outdir}/$counter.temp") or die "Cannot open $options{outdir}/$counter.temp";
	my @intemp = <INTEMP>;
	print OUTTABLE $intemp[0];
	close INTEMP;
	system "rm $options{outdir}/$counter.temp";
}

#############################################################################
# Subroutine to parse pfamscan outputs
#############################################################################

sub run_pfamscan_parser($){
	my $counter = $_[0];								
	my @return;
	if ($options{quiet} eq "N"){ print "Parsing domains from on $inpfamscans[$counter]: ", $counter + 1, " of ", scalar @inpfamscans, "\n"};

	my @domains;
	$inpfamscans[$counter] =~ s/\s*$//;
	(my $name = $inpfamscans[$counter]) =~ s/\.\w*$//; 
	$name =~ s/^.*\///;	
	open (OUTFILE, ">$options{outdir}/$name\_domains.faa") or die "Cannot open $options{outdir}/$name\_domains.faa";
	open (INFILE, $inpfamscans[$counter]) or die "Cannot open $inpfamscans[$counter]";
	my $dom_count = 0;
	while (<INFILE>){
		next if (/^#/);
		my @line = split /\s+/;
		next if (scalar @line != 15);
		$dom_count++; 
		my $locus_id   = $line[0];
		my $start      = $line[1];
		my $end        = $line[2];
		my $pfam_id    = $line[5];
		my $pfam_descr = $line[6];
		print OUTFILE ">$name\_$dom_count $locus_id $start-$end $pfam_id $pfam_descr\n";
		print OUTFILE substr($faas{$name}{$locus_id}, $start - 1, $end - $start + 1), "\n"; 
		push @domains, $pfam_id;
	}
	close OUTFILE;
	close INFILE;

	# collates summary lines as temp files to avoid thread collisions

	open (OUTSUMMARY, ">$options{outdir}/$counter.temp") or die "Cannot open $options{outdir}/$counter.temp";
	print OUTSUMMARY "$name\t", scalar @domains, "\t", join ",", @domains, "\n";
	close OUTSUMMARY;
}


