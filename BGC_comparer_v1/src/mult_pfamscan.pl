#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0 - January 6, 2018 - formal input parameters, more generalizable input, multithreaded
#
# script runs pfamscan.pl on a list of input faa files
# assumes that Pfam HMMs are located in the same place as the pfamscan script

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "mult_pfamscan.pl

DEPENDANCIES: Perl \"Parallel::ForkManager\" module

USAGE:
-c	number of cores to use						DEFAULT: 1		e.g., perl mult_pfamscan.pl -i inlist -c 4
-d	output directory for pfamscan output files files		DEFAULT: BGC_pfamscans	e.g., perl mult_pfamscan.pl -i inlist -d pfamscan_output
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl mult_pfamscan.pl -i inlist
-p	location of the pfamscan script and pfam HMMs			REQUIRED		e.g., perl mult_pfamscan.pl -i inlist -p ~/Tools/pfamscan
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl mult_pfamscan.pl -i inlist -q Y

OUTPUT FILES:
	a folder containing the pfamscan output for each file in the list of input files
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"d=s"  => \$options{outdir},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"p=s"  => \$options{pfamdir},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n $usage" if (!$options{infile});
die "Directory containing pfamscan.pl and its HMMS is not specified:\n $usage" if (!$options{pfamdir});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outdir}){   $options{outdir} = "BGC_pfamscans"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_pfamscan.pl	Jonathan Klassen	v2.0	Jan 6, 2018

parameters used:
	input file = $options{infile}
	directory containing pfamscan.pl and pfam HMMS = $options{pfamdir}\/
	output directory for pfamscan output files = $options{outdir}\/
	number of cores = $options{cores}
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

#############################################################################
# Multithreaded pfamscan
#############################################################################

# check if output directory exists

if (-d $options{outdir}){ die "$options{outdir}\/ already exists, existing mult_pfamscan.pl\n";}

# else make new output directory

mkdir $options{outdir} or die "Cannot make $options{outdir}";

# multithreaded pfamscan

my $pm = new Parallel::ForkManager($options{cores});
sub run_pfamscan($);
for my $a (0..$#infiles){
	$pm->start and next;
	my @return = run_pfamscan($a);
		$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# Subroutine to run pfamscan 
#############################################################################

sub run_pfamscan($){
	my $counter = $_[0];								# cluster id
	my @return;
	if ($options{quiet} eq "N"){ print "Running pfamscan.pl on $infiles[$counter]: ", $counter + 1, " of ", scalar @infiles, "\n"};

	# derive output file names from input file

	(my $name = $infiles[$counter]) =~ s/\.\w+$//;
	$name =~ s/^.+\///;

	# run pfamscan.pl

	system "perl $options{pfamdir}/pfam_scan.pl -fasta $infiles[$counter] -dir $options{pfamdir} -outfile $options{outdir}/$name.pfamscan";
}

