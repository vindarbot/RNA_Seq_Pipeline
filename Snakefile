import os
import sys
import re
import glob

include: "config.py"




FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))

CONDITIONS = list(set(x.split("_")[0] for x in SAMPLES))



DIRS = ['Reference','Reference/star/','Mapping','Mapping/Out','Trimming','featureCounts']

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
		"Reference/reference.fasta"



FASTA_NAME = os.path.basename(GET_GENOME)
GTF_NAME = os.path.basename(GET_GTF)


rule get_reference_files:
	output:
		"Reference/reference.fasta"

	params:
		get_genome = GET_GENOME,
		get_gtf = GET_GTF,
		get_description = GET_DESCRIPTION,
		fasta_name = FASTA_NAME,
		gtf_name = GTF_NAME

	message: ''' --- downloading fasta and gtf files --- '''

	shell: ''' 
		wget {params.get_genome}; mv {params.fasta_name} Reference/reference.fasta
		wget {params.get_gtf}; mv {params.gtf_name} reference.gff
		awk '{{ sub(/'ChrM'/,"mitochondria"); sub(/'ChrC'/,"chloroplast"); sub(/'Chr'/,"");print}}' reference.gff > reference_clean.gff
		rm reference.gff
		gffread reference_clean.gff -T -o reference.gtf
		rm reference_clean.gff
		wget {params.get_description}
		'''


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


GENOME = "Reference/reference.fasta"
GTF = "reference.gtf"


rule indexation_genome:
	input:
		genome = GENOME,
		gtf = GTF,
		starref = 'Reference/star/'

	output:
		"Reference/star/chrName.txt"

	threads: 16

	message: ''' --- Indexation du génome de référence --- '''

	shell: ' STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {input.starref} --genomeFastaFiles {input.genome}'


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

	threads: 16

	shell: ' STAR --runThreadN {threads} --sjdbGTFfile {input.gtf} --sjdbOverhang '+str(READ_LENGHT-1)+' --genomeDir {input.starref} \
	--outFileNamePrefix Mapping/{wildcards.sample} --readFilesIn {input.r1} {input.r2} --outSAMtype BAM SortedByCoordinate; \
	mv Mapping/{wildcards.sample}*.bam {output}; mv Mapping/*out* Mapping/Out '


rule sort_bam:
	input:
		"Mapping/{sample}.bam"

	output:
		"Mapping/{sample}.sorted.bam"

	threads: 16

	shell: ''' samtools sort {input} -o {output} -@ {threads} '''


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

	threads: 16

	shell: ''' featureCounts -T {threads} -p -t exon -g gene_id -a {params.gtf} -o {output} {input.mapping} '''





