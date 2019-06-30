#!/usr/bin/env python
import os
import sys
import re
import glob

configfile: "config.yaml"

RAWDATA_DIR = os.getcwd()

FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

wildcards = glob_wildcards('Experience/{fq_files}')

extension = [filename.split('.',1)[1] for filename in wildcards.fq_files][0]


if config["design"]["paired"]:

	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))

else:

	SAMPLES = list(set([ x.rstrip(extension) for x in FILES]))
	

CONDITIONS = list(set(x.split("_")[0] for x in SAMPLES))

CONDITION_TO_SAMPLES = {}

for condition in CONDITIONS:
	CONDITION_TO_SAMPLES[condition] = [sample for sample in SAMPLES if sample.split("_")[0] == condition]

DIRS = ['Reference','Reference/star/','Mapping','Mapping/Out','Trimming','featureCounts','DEG','logs']

for path in DIRS:
	if not os.path.exists(path):
		os.mkdir(path)


rule xpDesign: 		# Création d'un fichier txt qui décrit simplement le design expérimental, ceci est nécessaire pour l'étape d'analyse des gènes différentiellement exprimés sous R
	output:
		"experimentalDesign.txt"

	priority: 100

	run:
		with open("experimentalDesign.txt","w") as xpDesign:
			xpDesign.write("batch,condition\n")

			for condition,samples in CONDITION_TO_SAMPLES.items():
				for sample in samples:
					xpDesign.write(sample+".sorted.bam,"+condition+"\n")



genome = config["ref_files"]["genome"]
gtf = config["ref_files"]["gtf"]
description = config["ref_files"]["description"]
transcriptome = config["ref_files"]["transcriptome"]

GENOME = "Reference/reference.fasta"
GTF = "Reference/reference.gtf"
TRANSCRIPTOME = "Reference/transcriptome.fasta"


rule get_reference_files:	# Règle qui récupère le génome de référence ainsi que le fichier
							# d'annotation des gènes d'une espèce donnée
	input:
		"experimentalDesign.txt"

	output:
		fasta = "Reference/reference.fasta",
		gtf = "Reference/reference.gtf",
		transcripto = "Reference/transcriptome.fasta",
		description = "Reference/description.txt"

	params:
		get_genome = genome,
		get_gtf = gtf,
		get_description = description,
		fasta_name = os.path.basename(genome),
		gtf_name = os.path.basename(gtf),
		get_transcripto = transcriptome,
		transcripto_name = os.path.basename(transcriptome),
		description_name = os.path.basename(description)

	priority: 95

	message: ''' --- downloading fasta and gtf files --- '''

	shell: ''' 
		wget {params.get_genome}; mv {params.fasta_name} {output.fasta}
		wget {params.get_transcripto}; mv {params.transcripto_name} {output.transcripto}
		wget {params.get_gtf}; mv {params.gtf_name} reference.gff
		awk '{{ sub(/'ChrM'/,"mitochondria"); sub(/'ChrC'/,"chloroplast"); sub(/'Chr'/,"");print}}' reference.gff > reference_clean.gff
		rm reference.gff
		gffread reference_clean.gff -T -o {output.gtf}
		rm reference_clean.gff
		wget {params.get_description}
		mv {params.description_name} {output.description}
		'''



		