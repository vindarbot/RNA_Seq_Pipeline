#!/usr/bin/env python

##################################################
#  Fichier qui récupère plusieurs variables		 #
#  afin de générer l'analyse (fichier Snakefile) #
#  									 			 #
#  									 			 #
##################################################


### General

PAIRED_END = False

READ_LENGHT = 	150

GET_GENOME	=	"https://ics.hutton.ac.uk/atRTD/RTD2/AtRTDv2_QUASI_19April2016.fa"

GET_GTF		=	"https://ics.hutton.ac.uk/atRTD/RTD2/AtRTDv2_QUASI_19April2016.gtf"

GET_DESCRIPTION = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz"

### Trimming

TRIMMING =		True

# path to adapters
ADAPTERS =		"/happy_bin/bbmap/resources/adapters.fa"

# minimum length of reads to keep after trimming
minlen   =		25 

ktrim    =		"r" 
k        =		22 
qtrim    =		"rl" 
trimq    =		10 
hdist    =		1


### Mapping




### DEG Analysis

pADJ	=		0.11
LFC		= 		1