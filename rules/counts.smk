#!/usr/bin/env python
configfile: "config.yaml"

if config["design"]["paired"]:
	rule featureCounts_PE:
		input:
			gtf = "Reference/reference.gtf",
			mapping = expand("Mapping/{sample}.sorted.bam", sample=SAMPLES),
			index = expand("Mapping/{sample}.sorted.bam.bai", sample=SAMPLES)

		output:
			"featureCounts/counts.txt"

		priority: 75

		threads: 16

		shell: ''' featureCounts -p -s 2 -T {threads} -t exon -g gene_id -a {input.gtf} -o {output} {input.mapping} '''

else:
	rule featureCounts_SE:
		input:
			gtf = "Reference/reference.gtf",
			mapping = expand("Mapping/{sample}.sorted.bam", sample=SAMPLES),
			index = expand("Mapping/{sample}.sorted.bam.bai", sample=SAMPLES)

		output:
			"featureCounts/counts.txt"

		priority: 75

		threads: 16

		shell: ''' featureCounts -T {threads} -t exon -g gene_id -a {params.gtf} -o {output} {input.mapping} '''	

rule RPKM:
	input:
		"featureCounts/counts.txt"

	output:
		"DEG/RPKM.txt"

	shell: ''' python3 scripts/RPKM.py featureCounts/counts.txt'''