#!/usr/bin/env python

##############
# 	Script qui permet de générer des diagrammes de Venn avec des tailles de cercles proportionnels aux nombre de gènes
#	python3 scripts/venn.py path/to/genes_condition_1 path/to/genes_condition_2
#
#
##############
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import sys
import os


genes_condition_1 = [i.rstrip() for i in open(sys.argv[1]).readlines()]
genes_condition_2 = [i.rstrip() for i in open(sys.argv[2]).readlines()]
name_condition_1 = os.path.dirname(sys.argv[1]).split("/")[-1].replace("_"," VS ")
name_condition_2 = os.path.dirname(sys.argv[2]).split("/")[-1].replace("_"," VS ")

a = set(genes_condition_1)
b = set(genes_condition_2)

for gene in a.intersection(b):
	print(gene)

venn2([set(genes_condition_1), set(genes_condition_2)],set_labels = ('DASG','DEG'))
plt.title('ColHS VS hon4HS\n')
plt.show()


# venn2([set(genes_condition_1), set(genes_condition_2)],set_labels = ('DASG','DEG'))
# plt.title('Plantes contrôle VS plantes contrôle en stress thermique\n')