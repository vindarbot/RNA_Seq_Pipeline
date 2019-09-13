#!/usr/bin/env python
import os
configfile: "config.yaml"

# Alignement des echantillons pour l'analyse differentielle a l'echelle des genes.

if not os.path.exists('logs/Mapping'):
	os.mkdir('logs/Mapping')

rule indexation_genome:		# Indexation du génome de référence 
	input:
		gtf = "Reference/reference.gtf",
		genome = 'Reference/reference.fasta',
		starref = 'Reference/star/'

	output:
		"Reference/star/chrName.txt"

	log:
		"logs/Mapping/index.log"

	priority: 85

	threads: 8

	message: ''' --- Indexation du génome de référence --- '''

	shell: ' STAR --runThreadN {threads} --runMode genomeGenerate \
	--genomeDir {input.starref} \
	--genomeFastaFiles {input.genome} \
	--sjdbGTFfile {input.gtf} \
	--genomeSAindexNbases 11 \
	>{log} 2>&1'


if config["design"]["paired"]:
	rule mapping_PE:		
		input:
			index = "Reference/star/chrName.txt",
			starref = 'Reference/star/',
			r1 = 'Trimming/{sample}_R1.trim.'+config["extension"],
			r2 = 'Trimming/{sample}_R2.trim.'+config["extension"]

		output:
			"Mapping/{sample}.bam"

		priority: 80

		message: ''' --- Alignement des lectures --- '''

		threads: 4

		shell: ' STAR --runThreadN {threads} \
			--genomeDir {input.starref} \
			--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r1} {input.r2} \
			--readFilesCommand "gunzip -c" \
			--outSAMtype BAM Unsorted; \
			mv Mapping/{wildcards.sample}Aligned.out.bam {output};\
			mv Mapping/{wildcards.sample}*out* logs/Mapping'		


else:
	rule mapping_SE:		
		input:
			index = "Reference/star/chrName.txt",
			starref = 'Reference/star/',
			r = 'Trimming/{sample}.trim.'+config["extension"]

		output:
			"Mapping/{sample}.bam"

		priority: 80

		message: ''' --- Alignement des lectures --- '''

		threads: 4

		shell: ' STAR --runThreadN {threads} --genomeDir {input.starref} --readFilesCommand "gunzip -c" \
			--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r} --outSAMtype BAM Unsorted; \
			mv Mapping/{wildcards.sample}*.bam {output}; mv Mapping/*out* logs/Mapping '



rule sort_bam:
	input:
		"Mapping/{sample}.bam"

	output:
		"Mapping/{sample}.sorted.bam"

	threads: 4

	shell: ''' samtools sort {input} -o {output} -@ {threads} '''



rule index_bam:
	input:
		"Mapping/{sample}.sorted.bam"

	output:
		"Mapping/{sample}.sorted.bam.bai"

	shell: ''' samtools index {input} > {output} '''

	