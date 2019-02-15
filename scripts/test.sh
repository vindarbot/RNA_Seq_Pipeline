
#!/bin/bash

Samples=($(ls Samples/*R1*))

NB_SAMPLES=${#Samples[@]}

for i in `seq 1 "$NB_SAMPLES"`
do
	R1=$(ls Samples/SA"$i"*R1*)
	R2=$(ls Samples/SA"$i"*R2*)

	
	bbduk.sh in1=$R1 in2=$R2 out1=Data/Trimming/$(basename "$R1" .fastq).trim.fastq out2=Data/Trimming/$(basename "$R2" .fastq).trim.fastq \
	ref=/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa \
	minlen=25 ktrim=l k=22 qtrim=rl trimq=20 hdist=1 tpe tbo

done