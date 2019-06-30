#!/usr/bin/env python
import os
configfile: "config.yaml"

if not os.path.exists('logs/Trimming'):
	os.mkdir('logs/Trimming')

if config["design"]["paired"]:
	rule trimming_PE: 		
		input:
			xpDesign = 'experimentalDesign.txt',
			adapters = config["ref_files"]["adapters"],
			r1 = 'Experience/{sample}_R1.'+config["extension"],
			r2 = 'Experience/{sample}_R2.'+config["extension"]

		output:
			r1 = 'Trimming/{sample}_R1.trim.fastq.gz',
			r2 = 'Trimming/{sample}_R2.trim.fastq.gz'

		log:
			"logs/Trimming/{sample}.log"

		priority: 90

		message: ''' --- Trimming  --- '''

		shell: ' bbduk.sh in1="{input.r1}" in2="{input.r2}" out1="{output.r1}" out2="{output.r2}" \
			ref="{input.adapters}" minlen=25 ktrim=r k=22 qtrim=rl trimq=20 hdist=1 tpe tbo ziplevel=7 >{log} 2>&1'



else:
	rule trimming_SE: 		
		input:
			xpDesign = 'experimentalDesign.txt',
			adapters = config["ref_files"]["adapters"],
			r = 'Experience/{sample}.'+config["extension"]

		output:
			r = 'Trimming/{sample}.trim.fastq.gz'

		log:
			"logs/Trimming/{sample}.log"

		priority: 90

		message: ''' --- Trimming  --- '''

		shell: ' bbduk.sh in="{input.r}" out="{output.r}" \
			ref="{input.adapters}" minlen=25 ktrim=r k=22 qtrim=rl trimq=20 hdist=1 tpe tbo ziplevel=7 >{log} 2>&1 '