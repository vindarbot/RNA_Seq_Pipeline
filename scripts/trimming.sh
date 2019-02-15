#!/bin/bash

mkdir -p Data/Trimming

### VARIABLES ###

ADAPTERS="/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa"

. scripts/initialisation.sh

##################

### FONCTIONS ###

bbduk() {

	R1=$(ls "$EXP"/"$SAMPLE"/"$SAMPLE"_"$REPLICATE"_R1*)
	R2=$(ls "$EXP"/"$SAMPLE"/"$SAMPLE"_"$REPLICATE"_R2*)

	NameR1=$(basename "$R1" | cut -d"." -f1)
	NameR2=$(basename "$R2" | cut -d"." -f1)

	OUT1=Data/Trimming/"$NameR1".trim.fastq
	OUT2=Data/Trimming/"$NameR2".trim.fastq
		
	bbduk.sh in1="$R1" in2="$R2" out1="$OUT1" out2="$OUT2" \
	ref="$ADAPTERS" minlen=25 ktrim=r k=22 qtrim=rl trimq=10 hdist=1 tpe tbo 

}


###############


for SAMPLE in $SAMPLES
	do
		for REPLICATE in `seq 1 $(ls "$EXP"/"$SAMPLE"/*R1* | wc -l)`
		do

			bbduk &

		done

	done


