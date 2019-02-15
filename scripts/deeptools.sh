#!/bin/bash

bamCoverage -b Alignement/alignement_SA1.sorted.bam --outFileFormat bedgraph -o alignement_SA1.bed --normalizeUsing RPKM -p 16 
