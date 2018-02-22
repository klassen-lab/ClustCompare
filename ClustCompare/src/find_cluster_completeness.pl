#!/usr/bin/perl
#
# Jonathan Klassen
# v1.0 - June 22/17
# v1.1 - June 27/17 - tweaked to take display_id (LOCUS) instead of accessions
# v2.0 - February 21/18 - rewritten multithreaded, formal parameter passing
#
# Determines if BGCs are complete by their overlapping with contig ends
# Queries each input genome to get genome lengths, converts to fna in dir specified by -p/
# For each cluster in cluster_ids.table, uses nucmer to align cluster to genome from which it came to find locations
# requires BioPerl

use strict;
use warnings;
use Bio::SeqIO;
use Parallel::ForkManager;
use Getopt::Long;

############################################################################
# Processes input arguments
############################################################################

my $usage = "find_cluster_completeness.pl

DEPENDANCIES: nucmer, Perl \"Parallel::ForkManager\" module

USAGE:
-c	number of cores to use							DEFAULT: 1				e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j ../Results/nodes.tsv -c 4
-h	displays this usage statement (also using --help)
-i	input table matching source gbk files to the BGCs annotated on them	REQUIRED				e.g., perl find_cluster_completeness -i files.list -j ../Results/nodes.tsv 
-j	input node table to add fragmentation data to				REQUIRED 				e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j nodes.tsv 
-n	output node table, modified from -j input to contain fragmentation	DEFAULT: nodes.tsv			e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j ../Results/nodes.tsv -n nodes.table
-o	output table describing the alignment of each BGC to their source seq	DEFAULT: BGC_fragmentation.tsv		e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j ../Results/nodes.tsv -o BGC_fragmentataion.table
-p	output folder for alignment outputs					DEFAULT: ../Data/BGC_fragmentation	e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j ../Results/nodes.tsv -p BGC_fragmentation
-q	run quietly, i.e., no STDOUT (Y or N)					DEFAULT: N				e.g., perl find_cluster_completeness -i ../Results/BGC_file_lookup.tsv -j ../Results/nodes.tsv -q Y

OUTPUT FILES:
	The table specified by -o contains the alignment of each BGC to the contig on which wat annotated and its fragmentation status, as judged by the BGC ending at a contig edge. If specified, the node table is updated with the BGC fragmentation status such that these data can be annotated on Cytoscape BGC similarity network.
";

# input arguements

my %options = ();
GetOptions(
	"c=i"  => \$options{cores},
	"h"    => \$options{help},
	"help" => \$options{help},
	"i=s"  => \$options{inlookup},
	"j=s"  => \$options{innodes},
	"n=s"  => \$options{outnodes},
	"o=s"  => \$options{outfile},
	"p=s"  => \$options{outdir},
	"q=s"  => \$options{quiet}
);

# display usage statement if called

die $usage if ($options{help});

# required arguements

die "Input table linking BGC folders to source genomes is not specified:\n, $usage" if (!$options{inlookup});
die "Input node table is not specified:\n, $usage" if (!$options{innodes});

# defaut arguments

unless ($options{cores}){    $options{cores} = 1};
unless ($options{outnodes}){ $options{outnodes} = "nodes.tsv"};
unless ($options{outfile}){  $options{outfile} = "BGC_fragmentation.tsv"};
unless ($options{outdir}){   $options{outdir} = "../Data/BGC_fragmentation"};
unless ($options{quiet}){    $options{quiet} = "N"};

# mystery input flags not allowed

die "Unrecognized command line arguments: @ARGV\n" if ($ARGV[0]);

# checks correct parameter formatting

die "Unrecognized command line arguements: -q = $options{quiet}\n$usage" unless ($options{quiet} eq "Y" or $options{quiet} eq "N");

# print parameters unless -q flag selected

print "-----------------------------------------------------------------------------
find_cluster_completeness.pl	Jonathan Klassen	v2.0	Feb 21, 2018

parameters used:
	input looup table linking clusters to genomes = $options{inlookup}
	input node table = $options{innodes}
	output node table = $options{outnodes}
	output fragmentation table = $options{outfile}
	output directory for data files = $options{outdir}
	number of cores = $options{cores}
	quiet = $options{quiet}
-----------------------------------------------------------------------------
" if ($options{quiet} eq "N");

############################################################################
# Process input genome tables
############################################################################

open (INGENOMES, $options{inlookup}) or die "Cannot open table linking genomes to their annotated BGCs\n$usage";
my @lookup_file;
while (<INGENOMES>){
	s/\s+$//;
	s/^cluster//;
	push @lookup_file, $_;
}

my %lookup_table;
my @lookup_header = split /\t/, $lookup_file[0];
shift @lookup_file;
foreach (@lookup_file){
	my @lookupline = split /\t/, $_;
	for my $a (1..$#lookupline){
		$lookup_table{$lookupline[0]}{$lookup_header[$a]} = $lookupline[$a];
	}
}
close INGENOMES;

# create list of genomes to analyze, based on Source_genome_file field in input lookup table

my %genome_names;
foreach (keys %lookup_table){
	my $genome_name = $lookup_table{$_}{Source_genome_file};
	$genome_names{$genome_name} = "y";
}	

############################################################################
# Process input node table
############################################################################

open (INNODES, $options{innodes}) or die "Cannot open node table listing BGCs to analyze\n$usage";
my @node_file;
while (<INNODES>){
	s/\s+$//;
	s/^cluster//;
	push @node_file, $_;
}

my %node_table;
my @node_header = split /\t/, $node_file[0];
shift @node_file;
foreach (@node_file){
	my @nodeline = split /\t/, $_;
	for my $a (1..$#nodeline){
		$node_table{$nodeline[0]}{$node_header[$a]} = $nodeline[$a];
	}
}
close INNODES;

############################################################################
# make hash of all contigs and their end positions
############################################################################

unless (-d $options{outdir}){ mkdir $options{outdir} or die "Cannot mkdir $options{outdir}\n$usage"; }

my %contig_lengths;
my $count = 0;
print "\nParsing genomes:\n\n";
my @genome_list = keys %genome_names;
foreach my $genome (@genome_list){
	$genome =~ s/\s*$//;
	(my $annotated_file = $genome) =~ s/^.*\///;
	$annotated_file =~ s/\..*$//;

	$count++;
	print "Parsing $annotated_file: $count of ", scalar @genome_list, " genomes\n";
	
	# use BioPerl to open gbk object

	my $seqio_obj = Bio::SeqIO->new(-file => $genome, -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # goes to next sequences in object
		my $display_id = $seq_obj->display_name;
		my $length = $seq_obj->length;
		$contig_lengths{$annotated_file}{$display_id} = $length;
	}

	# output file

	my $new_seqio_obj = Bio::SeqIO->new(-file => ">$options{outdir}/$annotated_file.fna", -format => 'fasta');

	# use BioPerl to open gbk object

	$seqio_obj = Bio::SeqIO->new(-file => $genome, -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # goes to next sequences in object, i.e., 1st one
		$new_seqio_obj->write_seq($seq_obj);
	}
}

############################################################################
# convert all input clusters into fna files
############################################################################

$count = 0;
print "\nParsing clusters:\n\n";
foreach (keys %lookup_table){

	my $cluster_file = $lookup_table{$_}{Cluster_file};
	$cluster_file =~ s/\s*$//;
	(my $name = $cluster_file) =~ s/\.\w+$//;;

	$count++;
	print "Parsing $cluster_file: $count of ", scalar keys %lookup_table, " clusters\n";
	
	# output file

	my $new_seqio_obj = Bio::SeqIO->new(-file => ">$options{outdir}/cluster$_.fna", -format => 'fasta');

	# use BioPerl to open gbk object

	my $seqio_obj = Bio::SeqIO->new(-file => $cluster_file, -format => "genbank");
	while (my $seq_obj = $seqio_obj->next_seq){ # goes to next sequences in object, i.e., 1st one
		$new_seqio_obj->write_seq($seq_obj);
	}
}

############################################################################
# align clusters to their genomes using nucmer
############################################################################


print "\nAligning clusters:\n\n";

# multithreaded cluster parsing

my $pm = new Parallel::ForkManager($options{cores});
sub cluster_aligner($);
foreach my $cluster (sort {$a <=> $b} keys %lookup_table){
	$pm->start and next;
	my @return = cluster_aligner($cluster);
	$pm->finish;
}
$pm->wait_all_children();
my $counter = 0;

############################################################################
# parse nucmer alignments
############################################################################

my %alignments;
foreach my $cluster (sort {$a <=> $b} keys %lookup_table){
	print "Parsing alignments: $cluster of ", scalar keys %lookup_table, " clusters\n";

	my $genome_name = $lookup_table{$cluster}{Source_genome_file};
	$genome_name =~ s/\s*$//;
	(my $annotated_file = $genome_name) =~ s/^.*\///;
	$annotated_file =~ s/\..*$//;

	open (INCOORDS, "$options{outdir}/cluster$cluster\_vs_$annotated_file.coords") or die "$options{outdir}/cluster$cluster\_vs_$annotated_file.coords";
	my @coords = <INCOORDS>;
	my @coord_line = split /\s+/, $coords[5];

	$alignments{$cluster}{genome_name} = $annotated_file;
	$alignments{$cluster}{ref_start}   = $coord_line[1];
	$alignments{$cluster}{ref_end}     = $coord_line[2];
	$alignments{$cluster}{contig_name} = $coord_line[12];

	close INCOORDS;
}

############################################################################
# create alignment summary table
############################################################################

# create output table, including fragmentation status

open (OUTFILE, ">$options{outfile}") or die "Cannot open $options{outfile}\n$usage";
print OUTFILE "Cluster_name\tRef_genome_name\tRef_genome_contig\tRef_start\tRef_end\t5_fragmented\t3_fragmented\n";
foreach my $cluster (sort {$a <=> $b} keys %lookup_table){
	my $five_frag;
	my $three_frag;

	print OUTFILE "cluster$cluster\t",
		      "$alignments{$cluster}{genome_name}\t",
		      "$alignments{$cluster}{contig_name}\t",
		      "$alignments{$cluster}{ref_start}\t",
		      "$alignments{$cluster}{ref_end}\t";

	# determine 5' fragmentation status
	if ($alignments{$cluster}{ref_start} == 1){
		print OUTFILE "Y\t";
		$five_frag = 1;
	}
	else {
		print OUTFILE "N\t";
		$five_frag = 0;
	}

	# determine 3' fragmentation status

	if ($alignments{$cluster}{ref_end} == $contig_lengths{$alignments{$cluster}{genome_name}}{$alignments{$cluster}{contig_name}}){
		print OUTFILE "Y\n";
		$three_frag = 1;
	}
	else {
		print OUTFILE "N\n";
		$three_frag = 0;
	}
	my $total_frag = $five_frag . $three_frag;
	$node_table{$cluster}{Fragmentation} = $total_frag;
}
close OUTFILE;

#############################################################################
# Output modified node table 
#############################################################################

open (OUTNODES, ">$options{outnodes}") or die "Cannot open node table listing BGCs to edit\n$usage";
push @node_header, "Fragmentation";
print OUTNODES join "\t", @node_header, "\n";
foreach my $cluster (sort {$a <=> $b} keys %node_table){
	print OUTNODES "cluster$cluster";
	for my $a (1..$#node_header){
		print OUTNODES "\t", $node_table{$cluster}{$node_header[$a]};
	}
	print OUTNODES "\n";
}

#############################################################################
# Subroutine to perform alignments 
#############################################################################

sub cluster_aligner($){
	my $counter = $_[0];								# cluster id
	my @return;

	my $genome_name = $lookup_table{$counter}{Source_genome_file};
	$genome_name =~ s/\s*$//;
	(my $annotated_file = $genome_name) =~ s/^.*\///;
	$annotated_file =~ s/\..*$//;

	print "\nAligning cluster$counter to $annotated_file: $counter of ", scalar keys %lookup_table, " clusters\n\n";
	
	# do nucmer alignment
	system "nucmer $options{outdir}/$annotated_file.fna $options{outdir}/cluster$counter.fna -p $options{outdir}/cluster$counter\_vs_$annotated_file";
	system "delta-filter -1 $options{outdir}/cluster$counter\_vs_$annotated_file.delta > $options{outdir}/cluster$counter\_vs_$annotated_file.filtered.delta";
	system "show-coords $options{outdir}/cluster$counter\_vs_$annotated_file.filtered.delta > $options{outdir}/cluster$counter\_vs_$annotated_file.coords";
	return ($counter); 
}




