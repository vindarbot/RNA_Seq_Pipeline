#!/bin/python3

import sys
import os

### 
# Script qui prend en entrée un fichier de sortie de rMATS pour un évenèment de splicing donné, et ressort
# les évènement différentiellement exprimés entre les deux conditions.
# 
# On filtre d'abord le fichier pour éliminés les évè!nements donc le nombre moyen de reads associés entre les échantillons
# est inférieur à 5
#
# Esnuite, on garde les évènements dont la valeur de FDR (False Discovery Rate) est inférieur à 0.05, 
# et dont la valeur absolue du IncLevel est supérieur à 0.1
#
# Enfin, on sélectionne uniquement les colonnes intéressantes à garder.



rMATs_file = open(sys.argv[1],"r")

file_name = os.path.basename(sys.argv[1]).rstrip(".txt")

columnsToKeep = [0,1,3,4,5,6,7,8,9,10,12,14,19,22]


with open("topSplicingEvents/"+file_name+".top.txt","w") as rMATSs_top_file:

	for ligne in rMATs_file.readlines():

		if ligne.startswith("ID"):

			for i in columnsToKeep:

				rMATSs_top_file.write(ligne.split()[i]+"\t")
			rMATSs_top_file.write("\n")



		else:
			nbTotalReads = sum(list(map(int, ligne.split()[12].split(","))))+sum(list(map(int, ligne.split()[14].split(","))))

			nbSamples = len(ligne.split()[12].split(","))+len(ligne.split()[14].split(","))

			avgNbReads = nbTotalReads / nbSamples

			if avgNbReads > 6 and float(ligne.split()[19]) < 0.05 and abs(float(ligne.split()[22])) > 0.1 :


				for i in columnsToKeep:

					rMATSs_top_file.write(ligne.split()[i]+"\t")

				rMATSs_top_file.write("\n")

rMATs_file.close()










