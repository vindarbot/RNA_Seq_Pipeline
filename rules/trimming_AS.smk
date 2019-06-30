#!/usr/bin/env python
import os
import sys
import re
import glob

configfile: "config.yaml"

CONDITION1 = config["design"]["condition_1"]
CONDITION2 = config["design"]["condition_2"]

DIR_COMPARAISON = "DAS/"+CONDITION1+"_VS_"+CONDITION2

DIRS = ['Reference/star_1stpass','Trimming_AS','DAS','pass1','genomeForPass2','pass2',"DAS",DIR_COMPARAISON,DIR_COMPARAISON+"/rMATS_output/",DIR_COMPARAISON+"/topSplicingEvents","logs/Trimming_AS"]

for path in DIRS:
	if not os.path.exists(path):
		os.mkdir(path)



if config["design"]["paired"]:
	rule trimming_AS_PE: 		
		input:
			r1 = 'Trimming/{sample}_R1.trim.fastq.gz', 
			r2 = 'Trimming/{sample}_R2.trim.fastq.gz'

		output:
			r1 = 'Trimming_AS/{sample}_R1.trim.fastq.gz',
			r2 = 'Trimming_AS/{sample}_R2.trim.fastq.gz'

		params:
			cut_trim = config["DASG"]["cut_trim"]

		log:
			"logs/Trimming_AS/{sample}.log"

		priority: 65

		message: ''' --- Second Trimming Step (for DASG)  --- '''

		shell: ' cutadapt -l {params.cut_trim} -m {params.cut_trim} -o {output.r1} -p {output.r2} {input.r1} {input.r2} -j 8 --pair-filter=any >{log} 2>&1'


else:
	rule trimming_AS_SE: 		
		input:
			'Trimming/{sample}.trim.fastq.gz'

		output:
			'Trimming_AS/{sample}.trim.fastq.gz'

		params:
			cut_trim = config["DASG"]["cut_trim"]

		log:
			"logs/Trimming_AS/{sample}.log"

		priority: 65

		message: ''' --- Second Trimming Step (for DASG) --- '''

		shell: ' cutadapt -l {params.cut_trim} -m {params.cut_trim} -o {output} {input} -j 8 >{log} 2>&1'





