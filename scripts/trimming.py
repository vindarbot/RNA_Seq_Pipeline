#!/usr/bin/env python

import sys
import os
import glob

##
# Variables
cwd = os.getcwd()
print(cwd)

# os.makedirs("Data/Trimming")

Samples = os.listdir(cwd+"/Samples")

file1 = [i for i in glob.glob("Samples/*R1*")]
file2 = [i for i in glob.glob("Samples/*R2*")]

for i in range(len(file1)):
	print("".join(os.path.basename(file1[i]).split(".")[:-1]))

for i in range(len(file1)):

	os.system("bbduk.sh in1="name1[i] in2=name2[i] out1=Data/Trimming/${name1%.*}.trim.fastq \
# 	out2=Data/Trimming/${name2%.*}.trim.fastq ref=/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa \
# 	minlen=25 ktrim=l k=22 qtrim=rl trimq=20 hdist=1 tpe tbo")



# 	os.system(bbduk.sh in1=name1[i] in2=name2[i] out1=Data/Trimming/${name1%.*}.trim.fastq \
# 	out2=Data/Trimming/${name2%.*}.trim.fastq ref=/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa \
# 	minlen=25 ktrim=l k=22 qtrim=rl trimq=20 hdist=1 tpe tbo)


# sys.stdout.write("un test\n")