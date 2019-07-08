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


if not os.path.exists("DAS/rMATS_output"):
	os.mkdir("DAS/rMATS_output")

if not os.path.exists('logs/rMATS'):
	os.mkdir('logs/rMATS')

###



if config["design"]["paired"]:
	rule rMATS_PE:
		input:
			gtf = "Reference/reference.gtf",
			starRef = "genomeForPass2/",
			b1 = "b1.txt",
			b2 = "b2.txt",
			dir = "DAS/rMATS_output"


		output:
			"DAS/rMATS_output/SE.MATS.JCEC.txt"

		params:
			read_length = config["DASG"]["cut_trim"],
			rmats_dir = config['dir']['rmats'],
			strand = config["design"]["lib_type"]

		log:
			"logs/rMATS/rMATS.log"

		priority: 50

		message: ''' --- Second Trimming Step (for DASG)  --- '''

		shell:'''
			python2.7 {params.rmats_dir} --b1 {input.b1} --b2 {input.b2} --readLength {params.read_length} -t paired --gtf {input.gtf} --bi {input.starRef} --od {input.dir} --libType {params.strand} >{log} 2>&1''' 

else:
	rule rMATS_SE:
		input:
			gtf = "Reference/reference.gtf",
			starRef = "genomeForPass2/",
			b1 = "b1.txt",
			b2 = "b2.txt",
			dir = "DAS/rMATS_output"


		output:
			"DAS/rMATS_output/SE.MATS.JCEC.txt"

		params:
			read_length = config["DASG"]["cut_trim"],
			rmats_dir = config['dir']['rmats'],
			strand = config["design"]["lib_type"]

		log:
			"logs/rMATS/rMATS.log"

		priority: 50

		message: ''' --- Second Trimming Step (for DASG)  --- '''

		shell:'''
			python2.7 {params.rmats_dir} --b1 {input.b1} --b2 {input.b2} --readLength {params.read_length} -t single --gtf {input.gtf} --bi {input.starRef} --od {input.dir} --libType {params.strand} >{log} 2>&1''' 

