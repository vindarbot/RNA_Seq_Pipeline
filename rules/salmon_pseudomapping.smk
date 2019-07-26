#!/usr/bin/env python
import os

configfile: 'config.yaml'

if not os.path.exists('logs/Salmon'):
	os.mkdir('logs/Salmon')

DIRS = ['Salmon','Salmon/index','Salmon/quants','DTU']

for path in DIRS:
	if not os.path.exists(path):
		os.mkdir(path)	


rule salmon_index:
	input:
		"Reference/transcriptome.fasta"
	output:
		"Salmon/index/header.json"

	params:
		dir = "Salmon/index"

	threads:
		8
	log:
		"logs/Salmon/index.log"

	message:
		" --- Salmon indexing --- "

	shell:
		" salmon index -p {threads} -t {input} -i {params.dir} >{log} 2>&1"


if config["design"]["paired"]:

	rule salmon_mapping_PE:
		input:
			r1 = 'Trimming/{sample}_R1.trim.fastq.gz',
			r2 = 'Trimming/{sample}_R2.trim.fastq.gz',
			index = 'Salmon/index',
			outindex = "Salmon/index/header.json"

		output:
			'Salmon/quants/{sample}/quant.sf'

		params:
			boots = 10,
			outdir = 'Salmon/quants/{sample}'

		log:
			"logs/Salmon/quants_{sample}.log"

		shell:
			" salmon quant -l A --index {input.index} -1 {input.r1} -2 {input.r2} -o {params.outdir} --validateMappings --numBootstraps {params.boots} --writeMappings={params.outdir}/out.sam >>{log} 2>&1"

else:
	rule salmon_mapping_SE:
		input:
			r = 'Trimming/{sample}.trim.fastq.gz',
			index = 'Salmon/index',
			outindex = "Salmon/index/header.json"

		output:
			'Salmon/quants/{sample}/quant.sf'

		params:
			boots = 10,
			outdir = 'Salmon/quants/{sample}'

		log:
			"logs/Salmon/quants_{sample}.log"

		shell:
			" salmon quant -l A --index {input.index} -r {input.r} -o {params.outdir} --validateMappings --numBootstraps {params.boots} --writeMappings={params.outdir}/out.sam >>{log} 2>&1"






	