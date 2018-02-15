# ClustCompare
Pipeline to compare secondary metabolite biosynthetic gene clusters (BGCs)
--------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------

ClustCompare

Copyright Jonathan Klassen, 2018
This program is distributed under a GNU open-access license. See the accompanying LISCENSE file for terms of this license.

For further information, contact Jonathan Klassen: jonathan.klassen@uconn.edu

--------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------

Table of Contents:
1. Algorithm purpose and overview
2. Pipeline overview
3. Dependancies 
4. mult_antiSMASH pipeline
5. cluster_pfam_BBH_comparison pipeline
6. Visualizing the output data using Cytoscape
  
--------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------------------------------
1. Algorithm purpose and overview
--------------------------------------------------------------------------------------------------------------------------------

ClustCompare comprises several piplines designed to compare secondary metabolite biosynthetic gene clusters (BGCs) to each other based on their orthologous domain content. (In principle any gene cluster could be analyzed using this protocol, not just BGCs). This pipeline takes a list of biosynthetic gene clusters as input and outputs a similarity score that compares each cluster to each other based on their shared domain content and the similarly of these domains to each other. This table can subsequently used in downstream applications, e.g., to make a similarity network using Cytoscape. 

Besides the core BGC comparison pipeline, scripts are also included to find gene clusters in input genomes and to generate metadata tables that complement the similarity score described above. Together, these scripts provide a complete pipeline to annotate and compare BGCs in any input genomes.

--------------------------------------------------------------------------------------------------------------------------------
2. Pipeline overview
--------------------------------------------------------------------------------------------------------------------------------
The overall pipeline is designed to produce a network of BGC similarities. There are two overarching pipeline scripts (antiSMASH_annotation.sh and cluster_pfam_BBH_comparison.sh) which each wrap several other scripts that perform each step of the BGC comparison pipeline. Each script create a table of node and edge attributes, respectively, where the nodes are each individual BGC and the edges are the similarities between each BGC. The overall pipeline is designed to use either annotated or unannotated genomes as inputs, which it then uses to predict BGCs as the basis for network construction. However, other data types can easily be accommodated by running each step of the pipeline separately using the appropriate input files. Most pipeline steps are multithreaded for computational efficiency.

By default, the pipeline has the following directory structure, which is assumed by antiSMASH_annotation.sh and cluster_pfam_BBH_comparison.sh:
```
BGC_comparer_v1/
	|----- src/ (containing all pipeline scripts)
	|----- Data/ (folder for all data used by the pipeline) 
	|	|-----genomes/ (containing all input genomes for the pipeline)
	|----- Results/ (folder for all results produced by the pipeline)
	|----- README (this README file)
```
Additional folders and a CHANGELOG are produced by the antiSMASH_annotation.sh and cluster_pfam_BBH_comparison.sh scripts to record and reorganize pipeline outputs. All annotations and intermediate files are moved to Data/ and all results tables are moved to Results/ by default. These can be customized using the appropriate options for each individual pipeline script.

--------------------------------------------------------------------------------------------------------------------------------
3. Dependancies
--------------------------------------------------------------------------------------------------------------------------------

The following dependancies need to be installed before running the ClustCompare pipeline. See the installation instructions provided by each link.

- BioPerl: http://bioperl.org/INSTALL.html
- Perl "Parallel::ForkManager" module: http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm
- antiSMASH: http://docs.antismash.secondarymetabolites.org/install/
- pfamscan: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/ Note that the filepath leading to these files needs to be specified for mult_pfamscan.pl (and therefore also in cluster_pfam_BBH_comparison.sh).
- BLAST+: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

--------------------------------------------------------------------------------------------------------------------------------
4. mult_antiSMASH pipeline
--------------------------------------------------------------------------------------------------------------------------------

Overview: These scripts run annotate BGCs on a list of genomes specified by the user. Each script can either be run individually to accommodate preexisting analyses, or as a unit to run the entire pipeline on a set of genomes for which BGCs have not yet been annotated.

-----------------------
antiSMASH_annotation.sh

Wrapper script that annotates BGCs in a series of input genomes and calculates various metadata relating to these clusters for subsequent steps. Includes three steps: (1) mult_antiSMASH.pl to run antiSMASH on multiple genomes to annotate the BGCs that they contain; (2) cluster_renamer.pl to collate the clusters annotated by antiSMASH and rename them as unique objects for subsequent analyses; also tabulates metadata for each annotated cluster (i.e., the nodes in a final cluster similarity network produced by subsequent steps); (3) clusters_per_genome to count the number of clusters of each type present in each genome. antiSMASH_annotation.sh runs without parameters, looking for input genomes in ../Data/genomes (i.e., assumes script is run from src/ as downloaded). Accepts .gbk, .gbff, .embl, .fna, .fa, .fasta nucleotide files as unput. Will gunzip if required. Runs mult_antiSMASH using all cores. Move all annotations to ../Data/antiSMASH_annotations

---------------------
(1) mult_antiSMASH.pl

Takes as input a list of genomes to be annotated. This pipeline will run most efficiently using annotated input files (e.g., gbk or embl) and preserve any metadata associated with those input files, e.g., species names and descriptions. However, antiSMASH will annotate raw fa, fna, or fasta files using default metadata (albeit with a longer run time). By default does not include antiSMASH's knownclusterblast functionality to shorten runtime.
```
DEPENDANCIES: run_antiSMASH.py, Perl \"Parallel::ForkManager\" module

USAGE:
-c	number of cores to use					DEFAULT: 1		e.g., perl mult_antiSMASH.pl -i inlist -c 4
-h	displays this usage statement (also using --help)
-i	input list of paths to the gbk files to be analyzed	REQUIRED		e.g., perl mult_antiSMASH.pl -i inlist
-k	enable the antiSMASH knownclusterblast function		DEFAULT: N		e.g., perl mult_antiSMASH.pl -i inlist -k Y
-o	output file listing intput & output files		DEFAULT: files.list	e.g., perl mult_antiSMASH.pl -i inlist -o outputs.list
-q	run quietly, i.e., no STDOUT (Y or N)			DEFAULT: N		e.g., perl multantiSMASH.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder that contains the antiSMASH annotation for that genome. Folders are linked to their corresponding input file in the table specified by -o.```
```
----------------------
(2) cluster_renamer.pl

Takes as input a list of gbk files, each containing a different BGC as annotated by antiSMASH. In principle, will run on non-antiSMASH clusters too but will not be able to extract data relateding to cluster type unless this information is included as a "/product=" annotation of the global "cluster" primary annotation object. Renames all clusters so that they have non-overlapping names, which is required for each cluster to be included as a unique analysis object. Produces a metadata table that links each renamed cluster to their source data, source genome (identified by unique sequence DESCRIPTION annotations), and cluster type. This table contains the metadata that relates to the nodes of the cluster similarity networks produced by this pipeline.
```
DEPENDANCIES: Perl "Parallel::ForkManager" module, BioPerl

USAGE:
-c	number of cores to use						DEFAULT: 1		e.g., perl cluster_renamer.pl -i inlist -c 4
-d	output directory for renamed cluster gbk files			DEFAULT: BGCs		e.g., perl cluster_renamer.pl -i inlist -d renamed_clusters
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl cluster_renamer.pl -i inlist
-o	output tab-delimited table of cluster names (node table)	DEFAULT: nodes.tsv	e.g., perl cluster_renamer.pl -i inlist -o clusters.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder using the same name that contains the antiSMASH annotation for that genome
```
--------------------------
(3) clusters_per_genome.pl

Parses the node table produced by cluster_renamer.pl to produce a table counting the number of each cluster type in each genome. Running this script is not necessary for subsequent pipeline steps.
```
DEPENDANCIES: none

USAGE:
-h	displays this usage statement (also using --help)
-i	input node table produced by cluster_renamer.pl			REQUIRED				e.g., perl cluster_renamer.pl -i nodes.tsv
-o	output tab-delimited table of cluster types in each genome	DEFAULT: clusters_per_genome.tsv	e.g., perl cluster_renamer.pl -i nodes.tsv -o genome_clusters.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N				e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	for each input in the list specified by -i, a folder using the same name that contains the antiSMASH annotation for that genome
```
--------------------------------------------------------------------------------------------------------------------------------
4. cluster_pfam_BBH_comparison pipeline
--------------------------------------------------------------------------------------------------------------------------------

These scripts compare the orthologous pfam domain content between a set of input BGCs. The output of this pipeline is a table of BGC similarities, i.e., a table of edges in a cluster similarity networks that complements the node tables produced by the mult_antiSMASH pipeline described in section 2 above. 

-----------------------
cluster_pfam_BBH_comparison.sh

Wrapper script for the pipeline, consisting of 7 steps: (1) mult_gbk_to_faa.pl to convert gbk-formatted BGCs (e.g., as produced by antiSMASH) to multiple protein fasta files that can be searched by pfamscan; (2) mult_pfamscan.pl to perform a pfamscan analysis of domain content on each BGC; (3) mult_pfamscan_parser.pl to collate the domains present in each BGC and produces a protein fasta file for each domain; (4) mult_blastp.pl to compare each domain to all other domains in the entire dataset (including all BGCs); (5) blastp_to_BBH_list.pl to identify bidirection best BLAST hits (BBHs) for each domain as a means to identify orthologous domains; (6) make_cluster_similarities.pl to tabulate the number of orthologous domains that are shared between each cluster, out of all possible domains; and (7) cluster_similarity_table_parser.pl to filter the final edge table based on the number of shared domains, the percentage of shared domains out of all possible domains, and the sequence similarity between the shared domains. Designed to follow antiSMASH_annotation.sh, and assumes annotated BGCs in gbk format is located in ../Data/BGC_gbks/. All intermediate data is stored in ../Data/, and results tables are stored in ../Results/.

----------------------
(1) mult_gbk_to_faa.pl

Converts each gbk file listed in the input file list into a protein fasta (faa) file, the format needed for pfamscan.pl. Multithreaded.
```
DEPENDANCIES: Perl "Parallel::ForkManager" module, BioPerl

USAGE:
-c	number of cores to use						DEFAULT: 1		e.g., perl mult_gbk_to_faa.pl -i inlist -c 4
-d	output directory for BCG faa files				DEFAULT: BGC_faas	e.g., perl mult_gbk_to_faa.pl -i inlist -d renamed_faas
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl cluster_renamer.pl -i inlist
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl cluster_renamer.pl -i inlist -q Y

OUTPUT FILES:
	an output directory that contains protein faas for each gbk file in the input file list
```
----------------------
(2) mult_pfamscan.pl

Performs a pfamscan annotation of the domain content in each input faa file. Multithreaded. Requires explicitly specifying the folder in which the pfamscan script and its accompanying pfam HMMs (pressed using hmmpress) are located. 
```
DEPENDANCIES: Perl "Parallel::ForkManager" module, PfamScan and its accompanying HMMs executable and both in the same directory.

USAGE:
-c	number of cores to use						DEFAULT: 1		e.g., perl mult_pfamscan.pl -i inlist -c 4
-d	output directory for pfamscan output files files		DEFAULT: BGC_pfamscans	e.g., perl mult_pfamscan.pl -i inlist -d pfamscan_output
-h	displays this usage statement (also using --help)
-i	input list of paths to the cluster gbk files to be analyzed	REQUIRED		e.g., perl mult_pfamscan.pl -i inlist
-p	location of the pfamscan script and pfam HMMs			REQUIRED		e.g., perl mult_pfamscan.pl -i inlist -p ~/Tools/pfamscan
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N		e.g., perl mult_pfamscan.pl -i inlist -q Y

OUTPUT FILES:
	a folder containing the pfamscan output for each file in the list of input files
```
----------------------
(3) mult_pfamscan_parser.pl

Tabulates the number and identify of the pfam domains annotated in each BGC by pfamscan. Produces a separate multiple faa file for each BGC containing the sequence of each annotated domain in that BGC.
```
DEPENDANCIES: Perl "Parallel::ForkManager" module

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
```
----------------------
(4) mult_blastp.pl

Compares the pfam domains annotated in each BGC to all other domains annotated in all other BGCs included in the analysis using BLASTp. Requires a premade BLAST database containing all of the domains to be compared.
```
DEPENDANCIES: Perl "Parallel::ForkManager" module, blastp, premade BLAST databases of all BGC domains to be compared

USAGE:
-b	name of premade protein BLAST database			REQUIRED		e.g., perl mult_blastp.pl -i inlist -b all_domains.faa
-c	number of cores to use					DEFAULT: 1		e.g., perl mult_blastp.pl -i inlist -c 4
-d	output directory for blastp output files		DEFAULT: BGC_BLASTps	e.g., perl mult_blastp.pl -i inlist -b all_domains.faa -d blasts
-h	displays this usage statement (also using --help)
-i	input list of cluster domain faa files to be analyzed	REQUIRED		e.g., perl mult_blastp.pl -i inlist
-q	run quietly, i.e., no STDOUT (Y or N)			DEFAULT: N		e.g., perl mult_blastp.pl -i inlist -q Y

OUTPUT FILES:
	a folder containing the BLASTp results of the domains from each BGC compared to all domains
```
----------------------
(5) blastp_to_BBH_list.pl

Parses results of the BLASTp comparison of each domain to all other domains included in the analysis to identify BBH relationships between domains present in two different cluster, i.e., orthologous domains shared by these clusters. Note: default parameters tuned for relatively closely-related sequences (70% identical)
```
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
```
----------------------
(6) make_cluster_similarities.pl

Tabulate the number of domains that are shared between two BGCs, the proportion of domains in the clusters that are shared, and the percent identity shared by these domains. Includes all possible homologies between two clusters without filtering, to provide a master template on which cluster_similarity_parser.pl can be run using various settings.
```
DEPENDANCIES: none

USAGE:
-h	displays this usage statement (also using --help)
-i	input BBH table output by blastp_to_BBH.pl				REQUIRED				e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv
-j	input table of domains per cluster output by mult_pfamscan_parser.pl	REQUIRED				e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv
-o	output tab-delimited table of domain orthologs shared by two clusters	DEFAULT: raw_cluster_similarities.tsv	e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv -o BBHs.tsv
-q	run quietly, i.e., no STDOUT (Y or N)					DEFAULT: N				e.g., perl make_cluster_similarities.pl -i domain_BBHs.tsv -j cluster_domains.tsv -q Y

OUTPUT FILES:
	a table listing the number of domains that are shared between two BGCs (out of all possible domains) and how closely these domains are related to each other; this table contains all possible relationships without filtering
```
----------------------
(7) cluster_similarity_table_parser.pl

Filters the raw output table from make_cluster_similarities.pl based on: (i) the minimum cluster similarity score (0 low, 1 high); (ii) the minimum average ortholog % identity; (iii) the minimum number of shared domains; (iv) the minimum % of shared domains. Note that setting the minimum average ortholog % identity to less than the % identity threshold used in blastp_to_BBH_list.pl will not cause further filtering.
```
DEPENDANCIES: none

USAGE:
-a	minimum cluster similarity score threshold (range: 0-1)				DEFAULT: 0.3
-b	minimum average ortholog percent id threshold					DEFAULT: 70
-c	minimum number of shared domains between clusters				DEFAULT: 2
-d	minimum percent of shared domains between clusters				DEFAULT: 50

-h	displays this usage statement (also using --help)
-i	input cluster similariy table produced by make_cluster_similarities.pl		REQUIRED		e.g., perl clusters_similarly_table_parser.pl -i raw_edges.tsv
-o	output tab-delimited table of parsed cluster similarities (network edges)	DEFAULT: edges.tsv	e.g., perl clusters_similarly_table_parser.pl -i raw_edges.tsv -o cluster_edges.tsv
-q	run quietly, i.e., no STDOUT (Y or N)				DEFAULT: N				e.g., perl clusters_similarly_table_parser.pl -i inlist -q Y

OUTPUT FILES:
	a table listing the number of domains that are shared between two BGCs (out of all possible domains) and how closely these domains are related to each other; this table contains only relationships passing the specified threshold filters
```
--------------------------------------------------------------------------------------------------------------------------------
6. Visualizing the output data using Cytoscape
--------------------------------------------------------------------------------------------------------------------------------

The network produced by this pipeline can easily be visualized using standard network visualization. The section describes how to use Cytoscape  to achieve this. 
1. Install Cytoscape from (http://www.cytoscape.org/download.php).
2. Import the edge.tsv file as a new network.
	a. In the menu bar, select File/Import/Network/File
	b. Select "edges.tsv" or equivalent
	c. Click on the "Cluster_1" column header and change the "Meaning" tab to "Source Node"
	d. Click on the "Cluster_2" column header and change the "Meaning" tab to "Target Node"
	e. Click "OK" - this will produce a network linking each BGC to those that are similar to it. The attributes used to weight these edges can be changed using the "Layout" menu options.
3. Import the nodes.tsv file as network node attributes for this network
	a. In the menu bar, select File/Import/Table/file
	b. Select "nodes.tsv" or equivalent
	c. "Cluster_id" should be selected as the key by default, matching the node labels in edges.tsv
	d. Click "OK" - this will add the data in nodes.tsv as attributes that can be used to label the network, e.g., using the "Style" tab on the control panel.
