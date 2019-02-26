#!/usr/bin/env python

##################################################
#  Fichier qui récupère plusieurs variables		 #
#  afin de générer l'analyse (fichier Snakefile) #
#  									 			 #
#  									 			 #
##################################################


### General

PAIRED_END = False

READ_LENGHT = 	100

GET_GENOME	=	"https://ics.hutton.ac.uk/atRTD/RTD2/AtRTDv2_QUASI_19April2016.fa"

GET_GTF		=	"https://ics.hutton.ac.uk/atRTD/RTD2/AtRTDv2_QUASI_19April2016.gtf"

GET_DESCRIPTION = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz"

### Trimming

TRIMMING =		True


# rMATS demande en entrée une longeur de reads spécifique, on trim donc tous les reads de 13bp en 3' (même les reads n'ayant pas d'adaptapteurs)
CUT_TRIM   =		13

READ_LENGHT_TRIM = READ_LENGHT-CUT_TRIM



### Mapping




### DEG Analysis

pADJ	=		0.11
LFC		= 		1