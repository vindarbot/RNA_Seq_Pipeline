#!/usr/bin/env python
import os
import sys
import re
import glob

configfile: "config.yaml"

if config["design"]["paired"]:
	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))
else:
	SAMPLES = list(set([ x.rstrip(extension) for x in FILES]))

if not os.path.exists('logs/sign_DAS'):
	os.mkdir('logs/sign_DAS')

EVENTS = glob.glob("DAS/rMATS_output/*MATS.JCEC.txt")


rule rMATS_sig:
	input:
		"DAS/rMATS_output/SE.MATS.JCEC.txt"

	output:
		"DAS/DAS.txt"

	params:
		padj = config["DASG"]["pADJ"],
		psi = config["DASG"]["cutoff"]

	log:
		"logs/sign_DAS/sign_DAS.log"

	priority: 45

	shell:
		"python3 scripts/rMATS_filt.py -p {params.padj} -c {params.psi} >{log} 2>&1"