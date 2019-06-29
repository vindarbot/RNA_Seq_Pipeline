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

DIR_COMPARAISON = "DAS/"+config["design"]["condition_1"]+"_VS_"+config["design"]["condition_2"]

EVENTS = glob.glob(DIR_COMPARAISON+"/rMATS_output/*MATS.JCEC.txt")


rule rMATS_sig:
	input:
		DIR_COMPARAISON+"/rMATS_output/SE.MATS.JCEC.txt"

	output:
		DIR_COMPARAISON+"/topSplicingEvents/DAS.txt"

	params:
		padj = config["DASG"]["pADJ"],
		psi = config["DASG"]["cutoff"]

	log:
		"logs/sign_DAS/sign_DAS.log"

	run:
		for AS_EVENT in EVENTS:
			os.system("python3 scripts/rMATS_filt.py -p {params.padj} -c {params.psi} -e "+AS_EVENT+" >{log} 2>&1")