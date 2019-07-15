#!/usr/bin/env python

#############
#
#
#
#
#############

import sys
import os
import re
import subprocess
import glob
from collections import defaultdict
from decimal import *


file1 = open(sys.argv[1],"r")

file2 = open("DEG/RPKM.txt","r")

file3 = open("CHROMATINE-STATES/genes_to_states.txt",'r')


up = file1.readlines()

RPKM = file2.readlines()

STATES = file3.readlines()

RPKM_to_join = {}

for ligne in RPKM:

	RPKM_to_join[ligne.split()[0]] = ligne.split()[-1]

formate_states = {}

for ligne in STATES[1:]:

	formate_states[ligne.split('\t')[0].rstrip()] = ligne.split('\t')[1].rstrip()


all_states = []

for gene,states in formate_states.items():

	for state in states.split():

		if int(state) not in all_states:

			all_states.append(int(state))

			all_states.sort()


new_file = []

header= ["tair_locus","tair_symbol","log2FoldChange","padj","Description","mean_FPKM"]
for state  in all_states:
	header.append(str(state))


new_file.append(header)


for ligne in up[1:]:

	if ligne.startswith('tair_locus'):
		print(ligne)
		new_file.append(ligne.split('\t'))
	

	elif ligne.split()[0].rstrip() in RPKM_to_join:


		full_RPKM = ligne.rstrip().split('\t')
		full_RPKM.append(RPKM_to_join[ligne.split('\t')[0].rstrip()])

		if ligne.split()[0].rstrip() in formate_states:

			states_gene = [str(x) for x in formate_states[ligne.split()[0].rstrip()].split()]

			for state in header[6:]:
				if state in states_gene:
					full_RPKM.append('V')
				else:
					full_RPKM.append('X')



			# for state in formate_states[ligne.split()[0].rstrip()].split():
			# 	print(state)

		new_file.append(full_RPKM)


write_file = open(sys.argv[2],"w")

for ligne in new_file:
	write_file.write("\t".join(ligne))
	write_file.write("\n")


# for ligne in up:
# 	if ligne.split()[0].rstrip() in formate_states:
# 		i+=1
# 		print(i)



