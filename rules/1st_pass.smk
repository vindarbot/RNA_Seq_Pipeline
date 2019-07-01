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

if not os.path.exists('logs/1st_pass'):
	os.mkdir('logs/1st_pass')

###


rule genome_for_pass1:		# Indexation du génome de référence 
	input:
		genome = "Reference/reference.fasta",
		gtf = "Reference/reference.gtf",
		starref = 'Reference/star_1stpass/',

# 
	output:
		"Reference/star_1stpass/chrName.txt"

	params:
		read_length = config["DASG"]["cut_trim"] - 1

# 
	threads: 8

	priority: 60
# 
	message: ''' --- Indexation du génome de référence --- '''
# 
	shell: ' STAR --runThreadN {threads} --runMode genomeGenerate --genomeSAindexNbases 4 --sjdbOverhang {params.read_length} \
	--genomeDir {input.starref} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf}'



if config["design"]["paired"]:
	rule pass1_PE:
		input:
			gtf = "Reference/reference.gtf",
			genome = "Reference/reference.fasta",
			r1 = 'Trimming_AS/{sample}_R1.trim.'+config["extension"],
			r2 = 'Trimming_AS/{sample}_R2.trim.'+config["extension"],
			index_info = 'Reference/star_1stpass/chrName.txt'

		output:
			"pass1/{sample}SJ.out.tab"

		params:
			direct = "Reference/star"

		threads: 1

		shell:' STAR --runThreadN {threads} --genomeDir {params.direct} \
			--outFileNamePrefix pass1/{wildcards.sample} --readFilesIn {input.r1} {input.r2} \
			--outSAMtype BAM Unsorted --readFilesCommand "gunzip -c"'

else:
	rule pass1_SE:
		input:
			gtf = "Reference/reference.gtf",
			genome = "Reference/reference.fasta",
			r = 'Trimming_AS/{sample}.trim.'+config["extension"],
			index_info = 'Reference/star_1stpass/chrName.txt'

		output:
			"pass1/{sample}SJ.out.tab"

		params:
			direct = "Reference/star"

		threads: 4

		shell:' STAR --runThreadN {threads} --genomeDir {params.direct} \
			--outFileNamePrefix pass1/{wildcards.sample} --readFilesIn {input.r} \
			--readFilesCommand "gunzip -c" \
			--outSAMtype BAM Unsorted'




rule filt_SJ_out:
	input: 
		expand("pass1/{sample}SJ.out.tab", sample=SAMPLES)

	output:
		"SJ.filt.tab"

	shell:''' 
	gawk '$6==1 || ($6==0 && $7>2)' {input} >> {output};

	'''



rule SJ_db:
	input:
		"SJ.filt.tab"

	output: 
		"Reference/SJ.db"

	shell: '''
	awk 'BEGIN {{OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";}} {{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' \
	 {input} > {output};
	 mv {input} logs/1st_pass/
	'''