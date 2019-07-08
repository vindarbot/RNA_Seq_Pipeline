#!/usr/bin/env python
configfile: "config.yaml"

FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

if config["design"]["paired"]:
	SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))
else:
	SAMPLES = list(set([ x.rstrip(extension) for x in FILES]))

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

		shell: ''' featureCounts -T {threads} -t exon -g gene_id -a {input.gtf} -o {output} {input.mapping} '''	

rule FPKM:
	input:
		"featureCounts/counts.txt"

	output:
		"DEG/RPKM.txt"

	shell: ''' python3 scripts/RPKM.py featureCounts/counts.txt'''