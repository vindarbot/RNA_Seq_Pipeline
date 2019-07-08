#!/bin/python3

import sys
import os
import argparse
import glob
### 
# Script qui prend en entrée un fichier de sortie de rMATS pour un évenèment de splicing donné, et ressort
# les évènement différentiellement exprimés entre les deux conditions.
# 
# On filtre d'abord le fichier pour éliminés les évè!nements donc le nombre moyen de reads associés entre les échantillons
# est inférieur à 5
#
# Esnuite, on garde les évènements dont la valeur de FDR (False Discovery Rate) est inférieur à 0.05, 
# et dont la valeur absolue du IncLevelDifference est supérieur à 0.1
#
# Enfin, on sélectionne uniquement les colonnes intéressantes à garder.

parser = argparse.ArgumentParser(description='filter signficatives splice events')

parser.add_argument("-p","--padj",action='store',
	help="padj cutoff\n")

parser.add_argument("-c","--psi", action='store', 
	help="PSI cutoff")

args = parser.parse_args()

PSI = float(args.psi)
padj = args.padj


columnsToKeep = [0,1,3,4,5,6,7,8,9,10,12,13,14,15,19,22]

EVENTS = glob.glob("DAS/rMATS_output/*MATS.JCEC.txt")

for event in EVENTS:

	rMATs_file = open(event,"r")

	file_name = os.path.basename(event).rstrip(".txt")

	with open("DAS/topSplicingEvents/"+file_name+".top.txt","w") as rMATSs_top_file:

		for ligne in rMATs_file.readlines():

			if ligne.startswith("ID"):

				for i in columnsToKeep:

					rMATSs_top_file.write(ligne.split()[i]+"\t")
				rMATSs_top_file.write("\n")



			else:
				# On compte le nombre total de reads dans les 2 conditions (ligne.split()[12,13] pour accéeder aux nombre de reads des)
				# échantillons de la condition 1, et ligne.split()[14,15] pour accéder aux reads de la condition 2
				nbTotalReads = sum(list(map(int, ligne.split()[12].split(","))))+sum(list(map(int, ligne.split()[14].split(","))))+sum(list(map(int, ligne.split()[13].split(","))))+sum(list(map(int, ligne.split()[15].split(","))))

				# On compte le nombre d'échantillon , par exemple dans la colonne 12, si on a 3 échantillons on aura cette valeur la:
				# 12,14,8
				# Ainsi, en splittant la ligne par (",") et en comptant le nombre d'éléments obtenus, on a accès au nombre d'échantillons
				# On somme pour les 2 colonnes associés aux 2 conditiosns.
				nbSamples = len(ligne.split()[12].split(","))+len(ligne.split()[14].split(","))

				avgNbReads = nbTotalReads / nbSamples


				try:
					if avgNbReads > 10  and abs(float(ligne.split()[22])) > PSI :

						print(abs(float(ligne.split()[22])))
						for i in columnsToKeep:

							rMATSs_top_file.write(ligne.split()[i]+"\t")

						rMATSs_top_file.write("\n")
				except:
					break

	rMATs_file.close()


if len(os.listdir("DAS/topSplicingEvents")) == 5:

	liste_DAS = []

	for file in os.listdir("DAS/topSplicingEvents"):

		with open("DAS/topSplicingEvents/"+file,"r") as file:

			for ligne in file.readlines():

				if ligne.split()[1].startswith("\""):

					liste_DAS.append(ligne.split()[1].strip("\""))

	liste_unique_DAS = list(set(liste_DAS))


	with open("DAS/DAS.txt","w") as DAS:
		for gene in liste_unique_DAS:
			DAS.write(gene+"\n")









