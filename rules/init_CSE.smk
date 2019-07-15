#!/usr/bin/env python
import os

configfile: "config.yaml"

DIR = os.getcwd()


rule get_CSE:
	output:
		"CHROMATINE-STATES/README.md"
	shell: 
		'git clone https://github.com/Darbinator/CHROMATINE-STATES.git'

rule run_CSE:
	input:
		liste_deg = config["CSE"]["path_list_deg_name"],
		git = "CHROMATINE-STATES/README.md"

	output:
		"CSE_results/genes_to_states.txt"

	params:
		dir = "CSE_results",
		padj = config["CSE"]["pADJ"],
		dir_html = DIR+"/CHROMATINE-STATES/html/index.html"

	message: ''' -- Running Chromatine States Enrichment Analysis ---  '''

	shell:
		'python3 CHROMATINE-STATES/python/fisher_chrState.py -l {input.liste_deg} -o {params.dir} -p {params.padj}; \
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
