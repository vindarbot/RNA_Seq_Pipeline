#!/usr/bin/env python

import os
import sys
import re
import glob

configfile: "config.yaml"

### 

extension = config["extension"]

if config["trimming"]["exec"]:
	FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

else:
	FILES = [ os.path.basename(x) for x in glob.glob("Trimming/*") ]

	FILES = [ file.replace('.trim','') for file in FILES]

	for file in FILES:
		os.system('touch Experience/'+file)


if config["design"]["paired"]:
	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))
else:
	SAMPLES = list(set([ x.rstrip('.'+extension) for x in FILES]))



###

def get_input(wildcards):
	input_list = ["Reference/reference.fasta"]

	if config["design"]["paired"]:
		for fq in expand("Trimming/{sample}_R1.trim."+extension, sample=SAMPLES):
			input_list.append(fq)
	else:
		for fq in expand("Trimming/{sample}.trim."+extension, sample=SAMPLES):
			input_list.append(fq)

	if config["DEG"]["exec"]:
		input_list.append("DEG/tair_ids.txt")

	if config["DTU"]["exec"]:
		input_list.append("DTU/DTU.txt")
		for quant in expand('Salmon/quants/{sample}/quant.sf', sample=SAMPLES):
			input_list.append(quant)

	if config["DASG"]["exec"]:

		if config["design"]["paired"]:
			for fq_AS in expand("Trimming_AS/{sample}_R1.trim."+extension, sample=SAMPLES):
				input_list.append(fq_AS)

		else:
			for fq_AS in expand("Trimming_AS/{sample}.trim."+extension, sample=SAMPLES):
				input_list.append(fq_AS)

		input_list.append("DAS/DAS.txt")

	if config["CSE"]["exec"]:
		input_list.append("CSE_results/CSE_DEG_UP.txt")

	return input_list


include: "rules/get_packages.smk"
include: "rules/get_ref_files.smk"

if config["trimming"]["exec"]:
	include: "rules/trimming.smk"

if config["DEG"]["exec"]:
	include: "rules/classic_mapping.smk"
	include: "rules/counts.smk"
	include: "rules/run_DESeq2.smk"



if config["DASG"]["exec"]:
	include: "rules/trimming_AS.smk"
	include: "rules/1st_pass.smk"
	include: "rules/2sd_pass.smk"
	include: "rules/run_rMATS.smk"
	include: "rules/sign_DAS.smk"

if config["DTU"]["exec"]:
	include: "rules/salmon_pseudomapping.smk"
	include: "rules/run_RATs.smk"

if config["CSE"]["exec"]:
	include: "rules/init_CSE.smk"
	

rule all:	
	input:
		get_input


		

		