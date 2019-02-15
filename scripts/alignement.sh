#!/bin/bash

### VARIABLES ###


ALIGNEUR="$1"

INDEX="GenReference/tair10"
GEN_REF=$(ls GenReference/TAIR10.fasta)
#################


. scripts/initialisation.sh


### FONCTIONS ###

Reference() {

	mkdir -p GenReference

	if [ ! -f "GenReference/TAIR10.fasta" ]; then

	wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

	# Pour obtenir la description des gènes d'arabidopsis

	wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/gene_description_20131231.txt.gz 

	gunzip gene_description_20131231.txt.gz 


	mv TAIR10_chr_all.fas GenReference/TAIR10.fasta


	fi
}


Indexation() {

	

	case $ALIGNEUR in

	hisat) hisat2-build -p 16 "$GEN_REF" GenReference/tair10	
		;;

	tophat) bowtie2-build --threads 16 -f "$GEN_REF" GenReference/tair10
		;;

	star) STAR --runThreadN 16 --runMode genomeGenerate --genomeDir GenReference/ \
--genomeFastaFiles "$GEN_REF"
		;;

	esac
}


Alignement() {

	mkdir -p Alignement/Out

	R1=$(ls Data/Trimming/"$SAMPLE"_"$REPLICATE"_R1*)
	R2=$(ls Data/Trimming/"$SAMPLE"_"$REPLICATE"_R2*)

	OUTSAM="Alignement/${SAMPLE}_${REPLICATE}.sam"
	OUTBAM="Alignement/${SAMPLE}_${REPLICATE}.bam"
	
	case $ALIGNEUR in

	hisat) hisat2 -x "$INDEX" -p 2 -1 $R1 -2 $R2 -S $OUTSAM
	;;

	tophat) tophat2 -p 2 -o Alignement/"$SAMPLE"_"$REPLICATE" "$INDEX" $R1 $R2 

	mv Alignement/"$SAMPLE"_"$REPLICATE"/accepted_hits.bam $OUTBAM
	;;

	star) STAR --runThreadN 2 --genomeDir GenReference/ --outFileNamePrefix Alignement/"$SAMPLE"_"$REPLICATE" --readFilesIn $R1 $R2 --outSAMtype BAM SortedByCoordinate
	mv Alignement/"$SAMPLE"_"$REPLICATE"*.bam $OUTBAM
	mv Alignement/*out* Alignement/Out
	;;

	bbmap) bbmap.sh in1="$R1" in2="$R2" out="$OUTSAM" ref="$GEN_REF"
	;;

	esac 

}


Samtools() {

	R1=$(ls Data/Trimming/"$SAMPLE"_"$REPLICATE"_R1*)
	R2=$(ls Data/Trimming/"$SAMPLE"_"$REPLICATE"_R2*)

	if [ ! -f "Alignement/${SAMPLE}_${REPLICATE}.bam" ]; then

		samtools view -b -S Alignement/"$SAMPLE"_"$REPLICATE".sam | \
	 	samtools sort - -o Alignement/"$SAMPLE"_"$REPLICATE".sorted.bam 

	 	samtools index Alignement/"$SAMPLE"_"$REPLICATE".sorted.bam

#	 	rm Alignement/*sam



	 else
	 	samtools sort "$SAMPLE"_"$REPLICATE".bam -o Alignement/"$SAMPLE"_"$REPLICATE".sorted.bam

	 	samtools index Alignement/"$SAMPLE"_"$REPLICATE".sorted.bam

	 fi


}


mkdir featureCounts

featureCounts -p -t exon -g gene_id -a TAIR10.gtf -o featureCounts/counts.txt Alignement/*bam

#################



# Indexation

# for i in `seq 1 "$NB_SAMPLES"`
# do

# 	Alignement &

# done

# wait

for SAMPLE in $SAMPLES
	do
		for REPLICATE in `seq 1 $(ls "$EXP"/"$SAMPLE"/*R1* | wc -l)`
		do

			Alignement 

		done

	done








## Commande pour modifier le fichier d'annotation avant que les noms de séquences correspondent à celle du génome de référence.

# awk '{ sub(/'Chr'/,""); sub(/'Mt'/,"mitochondria"); sub(/'Pt'/,"chloroplast");print}' TAIR10_modified_crwn_4.gtf > TAIR10.gtf



#nohup bamCoverage -b Alignement/alignement_SA7.sorted.bam -of bedgraph --normalizeUsing RPKM -p 16 -o alignement_SA7.bedGraph &