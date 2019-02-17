import os
import sys
import re
import glob

include: "config.py"




FILES = [ os.path.basename(x) for x in glob.glob("Experience/*") ] 

SAMPLES = list(set([ "_".join(x.split("_")[:2]) for x in FILES]))

CONDITIONS = list(set(x.split("_")[0] for x in SAMPLES))

ADAPTERS="/Users/vindarbo/happy_bin/bbmap/resources/adapters.fa"





rule all:
	input: expand("Trimming/{sample}_trim_R1.fq", sample=SAMPLES)

if TRIMMING:

	os.mkdir("Trimming")

	rule trimming:
	    input:
	    	adapters = ADAPTERS,
	        r1 = 'Experience/{sample}_R1.fq.gz',
	        r2 = 'Experience/{sample}_R2.fq.gz'

	    output:
	        r1 = 'Trimming/{sample}_trim_R1.fq',
	        r2 = 'Trimming/{sample}_trim_R2.fq'

	    message: ''' --- Trimming des données --- '''

	    shell: ''' bbduk.sh in1="{input.r1}" in2="{input.r2}" out1="{output.r1}" out2="{output.r2}" \
	    ref="{input.adapters}" minlen=25 ktrim=r k=22 qtrim=rl trimq=10 hdist=1 tpe tbo '''
