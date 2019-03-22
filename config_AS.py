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

GET_GTF		=	"https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.gff"

GET_DESCRIPTION = "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz"

### Trimming

TRIMMING =		True


# rMATS demande en entrée une longeur de reads spécifique, on supprime donc 36bp en 3' pour tous les reads (même les reads n'ayant pas d'adaptapteurs)
# 151-36 = 115pb (longueur de reads à indiquer à cutadapt)
CUT_TRIM   =		115



### Mapping



