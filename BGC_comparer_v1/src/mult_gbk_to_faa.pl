#!/usr/bin/perl
#
# Jonathan Klassen
# v2.0 - January 6, 2018 - formal input parameters, more generalizable input, multithreaded
#
# script to convert gbks to faa
# needs to have BioPerl installed

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;
use Bio::SeqIO;

############################################################################
# Processes input arguments
############################################################################

my $usage = "mult_gbk_to_faa.pl

DEPENDANCIES: Perl \"Parallel::ForkManager\" module, BioPerl

USAGE:
-c	number of cores to use						DEFAULT: 1		e.g., perl mult_gbk_to_faa.pl -i inlist -c 4
-d	output directory for BCG faa files				DEFAULT: BGC_faas	e.g., perl mult_gbk_to_faa.pl -i inlist -d renamed_faas
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl cluster_renamer.pl -i inlist
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	an output directory that contains protein faas for each gbk file in the input file list
";

# input arguements

my %options = ();
GetOptions(
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

die "Input file is not specified:\n, $usage" if (!$options{infile});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outdir}){   $options{outdir} = "BGC_faas"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
mult_gbk_to_faa.pl	Jonathan Klassen	v2.0	Jan 6, 2018

parameters used:
	input file = $options{infile}
	output directory for renamed BGC faas = $options{outdir}\/
	number of cores = $options{cores}
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

if (-d $options{outdir}){ die "$options{outdir}\/ already exists, existing mult_gbk_to_faa.pl\n";}

# else make new output directory

mkdir $options{outdir} or die "Cannot make $options{outdir}";

# multithreaded cluster parsing

my $pm = new Parallel::ForkManager($options{cores});
sub cluster_parser($);
for my $a (0..$#clusters){
	$pm->start and next;
	my @return = cluster_parser($a);
	$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

#############################################################################
# subroutine to change gbks to faa
#############################################################################

sub cluster_parser($){
	my $counter = $_[0];	
	print "Converting $clusters[$counter] to faa format: ", $counter + 1, " of ", scalar @clusters, "\n";
	my @return;

	# derive outfile name from input file
	
	(my $name = $clusters[$counter]) =~ s/\.\w+$//;
	$name =~ s/^.+\///;
	open (OUTFILE, ">$options{outdir}/$name.faa") or die "Cannot open $options{outdir}/$name.faa";

	# use BioPerl to open gbk object

	my $seqio_obj = Bio::SeqIO->new(-file => $clusters[$counter], -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # goes to next contig
		for my $feat_obj ($seq_obj->get_SeqFeatures){ # goes to next feature on contig
			my $primary_tag = $feat_obj->primary_tag;
			next unless ($primary_tag eq "CDS"); # only use CDS info
			for my $tag ($feat_obj->get_all_tags){ # goes through CDS features
				if ($tag eq "locus_tag"){
				        print OUTFILE ">", $feat_obj->get_tag_values($tag), "\n";
				}
				elsif ($tag eq "translation"){
		        		print OUTFILE $feat_obj->get_tag_values($tag), "\n";
				}
			}	
		}	
	}
}





