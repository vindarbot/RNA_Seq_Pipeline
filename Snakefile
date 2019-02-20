import os
import sys
import re
import glob

include: "config.py"




FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))

CONDITIONS = list(set(x.split("_")[0] for x in SAMPLES))



DIRS = ['Reference/star/','Mapping','Mapping/Out','Trimming','featureCounts']

for path in DIRS:
	if not os.path.exists(path):
		os.mkdir(path)



rule all:
	input:
		expand('Trimming/{sample}_R1.trim.fastq', sample=SAMPLES),
		expand("Mapping/{sample}.bam", sample=SAMPLES),
		expand("Mapping/{sample}.sorted.bam", sample=SAMPLES),
		expand("Mapping/{sample}.sorted.bam.bai", sample=SAMPLES),
		"featureCounts/counts.txt"



rule trimming:
	input:
	    adapters = ADAPTERS,
	    r1 = 'Experience/{sample}_R1.fastq.gz',
	    r2 = 'Experience/{sample}_R2.fastq.gz'

	output:
	    r1 = 'Trimming/{sample}_R1.trim.fastq',
	    r2 = 'Trimming/{sample}_R2.trim.fastq'

	message: ''' --- Trimming  --- '''

	shell: ' bbduk.sh in1="{input.r1}" in2="{input.r2}" out1="{output.r1}" out2="{output.r2}" \
	    ref="{input.adapters}" minlen='+str(minlen)+' ktrim='+ktrim+' k='+str(k)+' qtrim='+qtrim+' trimq='+str(trimq)+' hdist='+str(hdist)+' tpe tbo '



rule indexation_genome:
	input:
		genome = GENOME,
		gtf = GTF,
		starref = 'Reference/star/'

	message: ''' --- Indexation du génome de référence --- '''

	shell: ' STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.genome}'


rule mapping_PE:
	input:
		gtf = GTF,
		index = "Reference/star/chrName.txt",
		starref = 'Reference/star/',
		r1 = 'Trimming/{sample}_R1.trim.fastq',
		r2 = 'Trimming/{sample}_R2.trim.fastq'

	output:
		"Mapping/{sample}.bam"

	message: ''' --- Alignement des lectures --- '''

	shell: ' STAR --runThreadN 2 --sjdbGTFfile {input.gtf} --sjdbOverhang '+str(READ_LENGHT-1)+'--genomeDir {input.starref} \
	--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM SortedByCoordinate; \
	mv Mapping/{wildcards.sample}*.bam {output}; mv Mapping/*out* Mapping/Out '


rule sort_bam:
	input:
		"Mapping/{sample}.bam"

	output:
		"Mapping/{sample}.sorted.bam"

	shell: ''' samtools sort {input} -o {output} '''


rule index_bam:
	input:
		"Mapping/{sample}.sorted.bam"

	output:
		"Mapping/{sample}.sorted.bam.bai"

	shell: ''' samtools index {input} > {output} '''


rule featureCounts:
	input:
		mapping = expand("Mapping/{sample}.sorted.bam", sample=SAMPLES),
		index = expand("Mapping/{sample}.sorted.bam.bai", sample=SAMPLES)

	output:
		"featureCounts/counts.txt"

	params:
		gtf = GTF

	shell: ''' featureCounts -p -t exon -g gene_id -a {params.gtf} -o {output} {input.mapping} '''













