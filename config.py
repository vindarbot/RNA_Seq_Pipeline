#!/usr/bin/env python

##################################################
#  Fichier qui récupère plusieurs variables		 #
#  afin de générer l'analyse (fichier Snakefile) #
#  									 			 #
#  									 			 #
##################################################


### Trimming
TRIMMING =		True
ADAPTERS =		"/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa"
minlen   =		25 
ktrim    =		"r" 
k        =		22 
qtrim    =		"rl" 
trimq    =		10 
hdist    =		1


### Mapping