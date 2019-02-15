#!/usr/bin/env python

import sys
import os
import re
import subprocess
import glob

### VARIABLES ###

samples = os.listdir("Data/Trimming")
aligneur = sys.argv[1]

#################

### FONCTIONS ###




#################


### GENOME DE REFERENCE 

if not os.path.isdir('Alignement/GenReference'):
	os.makedirs("Alignement/GenReference")

if not os.path.exists("Alignement/GenReference/TAIR10.fasta"):
	subprocess.run(''' 
		wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas; \

		wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz ; \

		gunzip gene_description_20131231.txt.gz; \

		mv TAIR10_chr_all.fas Alignement/GenReference/TAIR10.fasta; \

		gunzip Alignement/GenReference/TAIR10.fasta.gz;
		 ''', shell = True)


GenRef = "Alignement/GenReference/TAIR10.fasta"

###

### INDEXATION

if aligneur == "hisat":

	subprocess.run("hisat2-build -p 16 "+str(GenRef)+" Alignement/GenReference/tair10")

elif aligneur == "tophat":

	subprocess.run("bowtie2-build -p 16 "+str(GenRef)+" Alignement/GenReference/tair10")


INDEX = "Alignement/GenReference/tair10"


### ALIGNEMENT

for i in range(1,int(len(samples)/2)+1):
	
	R1 = glob.glob("Data/Trimming/SA"+str(i)+"*R1*")
	R2 = glob.glob("Data/Trimming/SA"+str(i)+"*R2*")


	if aligneur == "hisat":

		subprocess.run("hisat2 -x "+INDEX+" -p 2 -1 "+R1[0]+" \
	 -2 "+R2[0]+" -o Alignement/alignement_SA"+str(i)+".sam &", shell=True)


	elif aligneur == "tophat":

		subprocess.run("tophat2 -p 2 -o Alignement/alignement_SA"+str(i)+".sam "+INDEX+" "+R1[0]+" "+R2[0]+" &", shell=True)

	os.system("wait")

	subprocess.run(" samtools view -b -S alignement_SA"+str(i)+".sam | \
		samtools sort - -o Alignement/alignement_SA"+str(i)+".sorted.bam ; \
		samtools"

	 samtools index Alignement/alignement_SA"+str(i)".sorted.bam"










