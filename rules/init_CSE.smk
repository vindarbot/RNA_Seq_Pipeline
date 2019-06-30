#!/usr/bin/env python
import os

configfile: "config.yaml"

if not os.path.exists("CHROMATINE-STATES"):

	rule get_CSE:
		output:
			"CHROMATINE-STATES/README.md"
		shell: 
			'git clone https://github.com/Darbinator/CHROMATINE-STATES.git'

rule run_CSE:
	input:
		"DEG/tair_ids.txt"

	output:
		"CSE_results/genes_to_states.txt"

	params:
		dir = "CSE_results"

	message: ''' -- Running Chromatine States Enrichment Analysis ---  '''

	shell:
		'python3 CHROMATINE-STATES/python/fisher_chrState.py -l {input} -o {params.dir}'
