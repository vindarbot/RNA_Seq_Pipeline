import os
configfile: "config.yaml"


if not os.path.exists('logs/Mapping'):
	os.mkdir('logs/Mapping')

rule indexation_genome:		# Indexation du génome de référence 
	input:
		genome = 'Reference/reference.fasta',
		starref = 'Reference/star/'

	output:
		"Reference/star/chrName.txt"

	log:
		"logs/Mapping/index.log"

	threads: 16

	message: ''' --- Indexation du génome de référence --- '''

	shell: ' STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.genome} >{log} 2>&1'


if config["design"]["paired"]:
	rule mapping_PE:		
		input:
			gtf = "Reference/reference.gtf",
			index = "Reference/star/chrName.txt",
			starref = 'Reference/star/',
			r1 = 'Trimming/{sample}_R1.trim.fastq.gz',
			r2 = 'Trimming/{sample}_R2.trim.fastq.gz'

		output:
			"Mapping/{sample}.bam"


		message: ''' --- Alignement des lectures --- '''

		threads: 4

		shell: ' STAR --runThreadN {threads} --sjdbGTFfile {input.gtf} \
			--genomeDir {input.starref} \
			--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r1} {input.r2} \
			--readFilesCommand "gunzip -c" \
			--outSAMtype BAM SortedByCoordinate; \
			mv Mapping/{wildcards.sample}Aligned.sortedByCoord.out.bam {output};\
			mv Mapping/{wildcards.sample}*out* logs/Mapping'		


else:
	rule mapping_SE:		
		input:
			gtf = "Reference/reference.gtf",
			index = "Reference/star/chrName.txt",
			starref = 'Reference/star/',
			r = 'Trimming/{sample}.trim.fastq.gz'

		output:
			"Mapping/{sample}.bam"

		message: ''' --- Alignement des lectures --- '''

		threads: 6

		shell: ' STAR --runThreadN {threads} --sjdbGTFfile {input.gtf} --genomeDir {input.starref} --readFilesCommand "gunzip -c" \
			--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r} --outSAMtype BAM SortedByCoordinate; \
			mv Mapping/{wildcards.sample}*.bam {output}; mv Mapping/*out* logs/Mapping '



rule sort_bam:
	input:
		"Mapping/{sample}.bam"

	output:
		"Mapping/{sample}.sorted.bam"

	threads: 6

	shell: ''' samtools sort {input} -o {output} -@ {threads} '''



rule index_bam:
	input:
		"Mapping/{sample}.sorted.bam"

	output:
		"Mapping/{sample}.sorted.bam.bai"

	shell: ''' samtools index {input} > {output} '''