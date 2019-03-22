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
# et dont la valeur absolue du IncLevelDifference est supérieur à 0.1
#
# Enfin, on sélectionne uniquement les colonnes intéressantes à garder.



rMATs_file = open(sys.argv[1],"r")

file_name = os.path.basename(sys.argv[1]).rstrip(".txt")

dirName = "/".join(os.path.dirname(sys.argv[1]).split("/")[:-1])

print(dirName)