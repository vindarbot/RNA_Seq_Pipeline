#!/usr/bin/env python

import os
import sys
import re
import glob

configfile: "config.yaml"

### 

FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

if config["design"]["paired"]:
	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))
else:
	SAMPLES = list(set([ x.rstrip(extension) for x in FILES]))

DIR_COMPARAISON = "DAS/"+config["design"]["condition_1"]+"_VS_"+config["design"]["condition_2"]

###


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
		bbmap = "scripts/BBMap/README.md",
		featureCounts = "scripts/subread-1.6.1/README.txt",
		deg = "DEG/tair_ids.txt",
		dtu = "DTU/DTU.txt",
		das = "DAS/ColHS_VS_HMGA/topSplicingEvents/DAS.txt",
		ref = "Reference/reference.fasta",
		cse = "CSE_results/genes_to_states.txt"


		

		