#!/usr/bin/perl
#
# Jonathan Klassen
# v2.1
# February 15, 2018
#
# script annotates all genomes in genome_list.list (e.g., generate using list_genomes.pl)
# unzips if needed, then zips again at end
# uses an abbreviated antiSMASH annotation that only produces gbffs for each cluster
# only accepts preannotated gbff files
# v1.1 removes extra files that antiSMASH makes that are not reused, e.g., website output, whole genomes
# v1.2 disables as many extra antiSMASH outputs as possible
# v2.0 formal input parameters, multithreaded
# v2.1 outputs list of input file names and their corresponding output directories

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
-c	number of cores to use					DEFAULT: 1		e.g., perl mult_antiSMASH.pl -i inlist -c 4
-h	displays this usage statement (also using --help)
-i	input list of paths to the gbk files to be analyzed	REQUIRED		e.g., perl mult_antiSMASH.pl -i inlist
-k	enable the antiSMASH knownclusterblast function		DEFAULT: N		e.g., perl mult_antiSMASH.pl -i inlist -k Y
-o	output file listing intput & output files		DEFAULT: files.list	e.g., perl mult_antiSMASH.pl -i inlist -o outputs.list
-q	run quietly, i.e., no STDOUT (Y or N)			DEFAULT: N		e.g., perl multantiSMASH.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder that contains the antiSMASH annotation for that genome. Folders are linked to their corresponding input file in the table specified by -o.
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"k=s"  => \$options{knownclusterblast},
	"o=s"  => \$options{outfile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file is not specified:\n, $usage" if (!$options{infile});

# defaut arguments

unless ($options{cores}){             $options{cores} = 1};
unless ($options{knownclusterblast}){ $options{knownclusterblast} = "N"};
unless ($options{outfile}){           $options{outfile} = "files.list"};
unless ($options{quiet}){             $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{knownclusterblast}\n$usage" unless ($options{knownclusterblast} eq "Y" or $options{knownclusterblast} eq "N");
die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_antiSMASH.pl	Jonathan Klassen	v2.1	Feb 15, 2018

parameters used:
	input file = $options{infile}
	knownclusterblast activated = $options{knownclusterblast}
	output file = $options{outfile}
	number of cores = $options{cores}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

#############################################################################
# Loads list of input files
#############################################################################

my @genomes;
my @names;
open (INLIST, $options{infile}) or die $usage;

while (<INLIST>){
	s/\s+$//;
	next unless ($_ =~ /\w/);
	push @genomes, $_;
	my $name = $_;
	$name =~ s/^.*\///;
	$name =~ s/\..*$//;
	push @names, $name;
}

#############################################################################
# Multithreaded running of antiSMASH
#############################################################################

open (OUTFILE, ">$options{outfile}") or die "Cannot open $options{outfile}";
print OUTFILE "Input_file\tantiSMASH_output_directory\n";

sub mult_antiSMASH($);

my $pm = new Parallel::ForkManager($options{cores});

for my $a (0..$#genomes){
	$pm->start and next;
	my @return = mult_antiSMASH($a);
	print OUTFILE "$genomes[$a]\t$names[$a]\n";
	$pm->finish;
}
$pm->wait_all_children();

#############################################################################
# Subroutine to run antiSMASH on each genome
#############################################################################

sub mult_antiSMASH($){
	my $counter = $_[0];
	if ($options{quiet} eq "N"){ print "Aligning $genomes[$counter]: ", $counter + 1, " of ", scalar @genomes, "\n"};
	if ($options{knownclusterblast} eq "Y"){
		system "run_antismash.py -c 1 --knownclusterblast --disable-embl --disable-svg --outputfolder $names[$counter] $genomes[$counter]";
	}
	else {
		system "run_antismash.py -c 1 --disable-embl --disable-svg --outputfolder $names[$counter] $genomes[$counter]";
	}
}
