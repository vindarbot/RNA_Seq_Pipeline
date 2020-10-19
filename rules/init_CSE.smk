#!/usr/bin/env python
import os

# Fichier qui permet d'abord de telecharger le projet ChromatnStateEnrichment via github
# Ensuite de lancer l'analyse a partir de la liste des gènes dérégulés lors de l'analyse differentielle a l'echelle des genes
# Enfin, creation d'un fichier de sortie qui resume  l'ensemble de l'analyse d'expression a l'echelle des genes.
configfile: "config.yaml"

DIR = os.getcwd()


rule get_CSE:
	output:
		"ChromatineStateEnrichment/README.md"
	shell: 
		'git clone https://github.com/Darbinator/ChromatineStateEnrichment.git'

rule run_CSE:
	input:
		liste_deg = config["CSE"]["path_list_deg_name"],
		git = "ChromatineStateEnrichment/README.md"

	output:
		"CSE_results/genes_to_states.txt"

	params:
		dir = "CSE_results",
		padj = config["CSE"]["pADJ"],
		dir_html = DIR+"/ChromatineStateEnrichment/html/index.html"

	message: ''' -- Running Chromatine States Enrichment Analysis ---  '''

	shell:
		'python3 ChromatineStateEnrichment/python/fisher_chrState.py -l {input.liste_deg} -o {params.dir} -p {params.padj}; \
		ln -s {params.dir_html} CSE_results/results.html'

rule format_DEG_CSE_results:
	input:
		up = 'DEG/genes_up.txt',
		down = 'DEG/genes_down.txt',
		FPKM = 'DEG/RPKM.txt',
		gene_to_states = "CSE_results/genes_to_states.txt"

	output:
		up = "CSE_results/CSE_DEG_UP.txt",
		down = "CSE_results/CSE_DEG_DOWN.txt"

	shell: '''python3 scripts/generate_results.py {input.up} {output.up};
		python3 scripts/generate_results.py {input.down} {output.down} '''
