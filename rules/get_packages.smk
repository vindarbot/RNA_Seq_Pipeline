#!/usr/bin/env python
import os
import sys
import re
import glob

configfile: "config.yaml"

rule get_featureCounts:
	output:
		"scripts/subread-1.6.1/README.txt"

	shell:
		"git clone https://github.com/torkian/subread-1.6.1.git; \
		mv subread-1.6.1/ scripts/"