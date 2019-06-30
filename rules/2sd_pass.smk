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

CONDITION1 = config["design"]["condition_1"]

DIR_COMPARAISON = "DAS/"+CONDITION1+"_VS_"+CONDITION2


if not os.path.exists('logs/2sn_pass'):
	os.mkdir('logs/2sn_pass')

###


rule genomeForPass2:
	input:
		genome = "Reference/reference.fasta",
		SJ = "Reference/SJ.db"

	output:
		"genomeForPass2/chrName.txt"

	params:
		genomedir = "genomeForPass2/",
		read_length = config["DASG"]["cut_trim"] - 1

	threads: 4

	priority: 55

	shell:' STAR --runThreadN {threads} --runMode genomeGenerate --sjdbOverhang {params.read_length} \
	--genomeDir {params.genomedir} --genomeFastaFiles {input.genome} --sjdbFileChrStartEnd {input.SJ}'




if config["design"]["paired"]:
	rule pass2_PE:
		input:
			genome = "Reference/reference.fasta",
			r1 = 'Trimming_AS/{sample}_R1.trim.'+extension,
			r2 = 'Trimming_AS/{sample}_R2.trim.'+extension,
			starRef = "genomeForPass2/chrName.txt"

		output:
			"pass2/{sample}.bam"

		params:
			read_length = config["DASG"]["cut_trim"] - 1,
			direct = "genomeForPass2/"


		threads: 4

		shell:' STAR --runThreadN {threads} --chimSegmentMin 2 --outFilterMismatchNmax 3\
	 		--alignIntronMax 299999 \
	 		--genomeDir {params.direct} \
	 		--outFileNamePrefix pass2/{wildcards.sample} --readFilesIn {input.r1} {input.r2} \
	 		--sjdbOverhang  {params.read_length} --readFilesCommand "gunzip -c" \
	 		--outSAMtype BAM SortedByCoordinate; \
	 		mv pass2/{wildcards.sample}Aligned.sortedByCoord.out.bam {output};'
# --readFilesCommand "gunzip -c" 
else:
	rule pass2_SE:
		input:
			genome = "Reference/reference.fasta",
			r = 'Trimming_AS/{sample}.trim.fastq.gz',
			starRef = "genomeForPass2/chrName.txt"

		output:
			"pass2/{sample}.bam"

		params:
			read_length = config["DASG"]["cut_trim"] - 1,
			direct = "genomeForPass2/"


		threads: 4

		shell:' STAR --runThreadN {threads} --chimSegmentMin 2 --outFilterMismatchNmax 3\
	 		--alignIntronMax 299999 \
	 		--genomeDir {params.direct} \
	 		--outFileNamePrefix pass2/{wildcards.sample} --readFilesIn {input.r} \
	 		--sjdbOverhang  {params.read_length} --readFilesCommand "gunzip -c" \
	 		--outSAMtype BAM SortedByCoordinate; \
	 		mv pass2/{wildcards.sample}Aligned.sortedByCoord.out.bam {output};'



rule make_b_files:
	input:
		expand("pass2/{sample}.bam", sample=SAMPLES)

	output:
		b1 = "b1.txt",
		b2 = "b2.txt"

	params:
		file1 = [cond1_bam for cond1_bam in glob.glob("pass2/*bam") if cond1_bam.startswith("pass2/"+CONDITION1)],
		file2 = [cond2_bam for cond2_bam in glob.glob("pass2/*bam") if cond2_bam.startswith("pass2/"+CONDITION2)]

		
		

	run:
		with open(output.b1,'w') as file1:
			file1.write(",".join(params.file1))

		with open(output.b2,'w') as file2:
			file2.write(",".join(params.file2))

