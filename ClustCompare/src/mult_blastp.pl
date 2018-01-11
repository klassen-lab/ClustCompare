#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0 - January 6, 2018 - formal input parameters, more generalizable input, multithreaded
#
# BLASTs domains from set of BGCs against a database of all domains (db created separately, batteries not included)
# Assumes BLASTp is in $PATH


use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "mult_blastp.pl

DEPENDANCIES: Perl \"Parallel::ForkManager\" module, blastp, premade BLAST databases of all BGC domains to be compared

USAGE:
-b	name of premade protein BLAST database			REQUIRED		e.g., perl mult_blastp.pl -i inlist -b all_domains.faa
-c	number of cores to use					DEFAULT: 1		e.g., perl mult_blastp.pl -i inlist -c 4
-d	output directory for blastp output files		DEFAULT: BGC_BLASTps	e.g., perl mult_blastp.pl -i inlist -b all_domains.faa -d blasts
-h	displays this usage statement (also using --help)
-i	input list of cluster domain faa files to be analyzed	REQUIRED		e.g., perl mult_blastp.pl -i inlist
-q	run quietly, i.e., no STDOUT (Y or N)			DEFAULT: N		e.g., perl mult_blastp.pl -i inlist -q Y

OUTPUT FILES:
	a folder containing the BLASTp results of the domains from each BGC compared to all domains
";

# input arguements

my %options = ();
GetOptions(
	"b=s"  => \$options{db},
	"c=i"  => \$options{cores},
	"d=s"  => \$options{outdir},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{infile},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input file list of domain faas is not specified:\n $usage" if (!$options{infile});
die "No BLASTp database specified:\n $usage" if (!$options{db});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outdir}){   $options{outdir} = "BGC_BLASTps"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_blastp.pl	Jonathan Klassen	v2.0	Jan 6, 2018

parameters used:
	input file = $options{infile}
	protein BLAST database = $options{db}
	output directory for BLASTp output files = $options{outdir}\/
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

if (-d $options{outdir}){ die "$options{outdir}\/ already exists, existing mult_blastp.pl\n";}

# else make new output directory

mkdir $options{outdir} or die "Cannot make $options{outdir}";

# multithreaded pfamscan

my $pm = new Parallel::ForkManager($options{cores});
sub run_blastp($);
for my $a (0..$#infiles){
	$pm->start and next;
	my @return = run_blastp($a);
		$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# Subroutine to run BLASTp 
#############################################################################

sub run_blastp($){
	my $counter = $_[0];
	my @return;
	if ($options{quiet} eq "N"){ print "Running BLASTp of $infiles[$counter] vs all domains: ", $counter + 1, " of ", scalar @infiles, "\n"};

	# derive output file names from input file

	(my $name = $infiles[$counter]) =~ s/\.\w+$//;
	$name =~ s/^.+\///;

	# run BLASTp

	system "blastp -query $infiles[$counter] -db $options{db} -evalue 1e-5 -out $options{outdir}/$name\_vs_$options{db}.blastp -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"";

}

