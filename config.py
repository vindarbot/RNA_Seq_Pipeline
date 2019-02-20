#!/usr/bin/env python

##################################################
#  Fichier qui récupère plusieurs variables		 #
#  afin de générer l'analyse (fichier Snakefile) #
#  									 			 #
#  									 			 #
##################################################


### General

READ_LENGHT = 	150

GET_GENOME	=	"https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas"

GET_GTF		=	"https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff"

GET_DESCRIPTION = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz"

### Trimming

TRIMMING =		True

# path to adapters
ADAPTERS =		"/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa"

# minimum length of reads to keep after trimming
minlen   =		25 

ktrim    =		"r" 
k        =		22 
qtrim    =		"rl" 
trimq    =		10 
hdist    =		1


### Mapping

# path to reference genome (default = TAIR10.fasta)
GENOME	=		'Reference/TAIR10.fasta'

# path to gtf file with transcripts (default = TAIR10.gtf)
GTF		=		'TAIR10.gtf'
