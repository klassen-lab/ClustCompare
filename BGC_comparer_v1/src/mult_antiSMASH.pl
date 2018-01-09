#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0
# December 23, 2017
#
# script annotates all genomes in genome_list.list (e.g., generate using list_genomes.pl)
# unzips if needed, then zips again at end
# uses an abbreviated antiSMASH annotation that only produces gbffs for each cluster
# only accepts preannotated gbff files
# v1.1 removes extra files that antiSMASH makes that are not reused, e.g., website output, whole genomes
# v1.2 disables as many extra antiSMASH outputs as possible
# v2.0 formal input parameters, multithreaded

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "mult_antiSMASH.pl

DEPENDANCIES: run_antiSMASH.py, Perl \"Parallel::ForkManager\" module

USAGE:
-c	number of cores to use					DEFAULT: 1	e.g., perl mult_antiSMASH.pl -i inlist -c 4
-h	displays this usage statement (also using --help)
-i	input list of paths to the gbk files to be analyzed	REQUIRED	e.g., perl mult_antiSMASH.pl -i inlist
-q	run quietly, i.e., no STDOUT (Y or N)			DEFAULT: N	e.g., perl multantiSMASH.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder using the same name that contains the antiSMASH annotation for that genome
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"q=s"   => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n, $usage" if (!$options{infile});

# defaut arguments

unless ($options{cores}){ $options{cores} = 1};
unless ($options{quiet}){ $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_antiSMASH.pl	Jonathan Klassen	v2.0	Dec 23, 2017

parameters used:
	input file = $options{infile}
	number of cores = $options{cores}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @genomes;
open (INLIST, $options{infile}) or die $usage;

while (<INLIST>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @genomes, $_;
}

#############################################################################
# Multithreaded running of antiSMASH
#############################################################################

sub mult_antiSMASH($);

my $pm = new Parallel::ForkManager($options{cores});

for my $a (0..$#genomes){
	$pm->start and next;
	my @return = mult_antiSMASH($a);
	$pm->finish;
}
$pm->wait_all_children();

#############################################################################
# Subroutine to run antiSMASH on each genome
#############################################################################

sub mult_antiSMASH($){
	my $counter = $_[0];
	if ($options{quiet} eq "N"){ print "Aligning $genomes[$counter]: ", $counter + 1, " of ", scalar @genomes, "\n"};
	system "run_antismash.py -c 1 --disable-embl --disable-svg $genomes[$counter]";
	
	# removes annotations files that are not needed

#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/css"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/css";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/images"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/images";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/js"){	
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/js";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/nrpspks_predictions_txt"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/nrpspks_predictions_txt";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/structures"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/structures";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/svg"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/svg";
#	}
#	if (-d "../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/txt"){
#		system "rm -r ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/txt";
#	}
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.xml";
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.final.*";
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.zip";
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.xls";
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.js";
#	system "rm ../Data/antiSMASH_annotations/raw_annotations/$genome_out_name/*.html";


	# zip input files, whether or not they started that way

#	system "gzip $genome_name";
}
