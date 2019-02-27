#!/usr/bin/env python

##################################################
#  Fichier qui récupère plusieurs variables		 #
#  afin de générer l'analyse (fichier Snakefile) #
#  									 			 #
#  									 			 #
##################################################


### General

RMATS = '/data/home/darbotv/happy_bin/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py'


PAIRED_END = True

READ_LENGHT = 	150

GET_GENOME	=	"https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas"

GET_GTF		=	"https://ics.hutton.ac.uk/atRTD/RTD2/AtRTDv2_QUASI_19April2016.gtf"

GET_DESCRIPTION = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz"

### Trimming

TRIMMING =		True


# rMATS demande en entrée une longeur de reads spécifique, on trim donc tous les reads de 13bp en 3' (même les reads n'ayant pas d'adaptapteurs)
CUT_TRIM   =		13

READ_LENGHT_TRIM = READ_LENGHT-CUT_TRIM

### Mapping



