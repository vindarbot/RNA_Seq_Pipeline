Installation
--------------
You may use Subread package by directly downloading a binary release suitable for your operating system (no compilation is needed), or you may build it from the source. Here we describe how to install it from the source.

Download the latest version of Subread package from http://subread.sourceforge.net/. The source release includes a keyword 'source' in the name of the tar ball. Uncompress the tar ball, enter the 'src' directory and issue the following command to build it for Linux OS :

make -f Makefile.Linux

For Mac OS, use command:

make -f Makefile.MacOS

For FreeBSD OS, use command:

gmake -f Makefile.FreeBSD

If the build is successful, a new directory called 'bin' will be created under the home directory of the package (ie. one level up from 'src' directory). The 'bin directory contains all the generated executables. To enable easy access to these executables, you may copy the executables to a system directory such as '/usr/bin' or add the path to the executables to your search path (add path to your environment variable `PATH').

Content
--------------
annotation    Directory including NCBI RefSeq gene annotations for genomes 'hg19', 'hg38', 'mm10' and 'mm9'.
              Each row is an exon. Entrez gene identifiers and chromosomal coordinates are provided for each exon.
bin           Directory including executables after compilation (or directly available from a binary release). 
doc           Directory including the users manual.
LICENSE       The license agreement for using this package.
README.txt    This file.
src           Directory including source code (binary releases do not have this directory).
test          Directory including test data and scripts.

A Quick Start
--------------
Build index for a reference genome:

  subread-buildindex -o my_index chr1.fa chr2.fa ...
  (You may provide a single FASTA file including all chromosomal sequences).

Align a single-end RNA-seq dataset to the reference genome:

  subread-align -i my_index -r reads.txt -t 0 -o subread_results.bam

Align a paired-end genomic DNA-seq dataset to the reference genome:

  subread-align -i my_index -r reads1.txt -R reads2.txt -t 1 -o subread_results_PE.bam

Detect exon-exon junctions from a paired-end RNA-seq dataset (read mapping results are also produced):

  subjunc -i my_index -r reads1.txt -R reads2.txt -o subjunc_results.bam

Assign mapped RNA-seq reads to mm10 genes using inbuilt annotation:

  featureCounts -a ../annotation/mm10_RefSeq_exon.txt -F 'SAF' -o counts.txt subread_results.bam

Assign mapped RNA-seq reads to hg38 genes using a public GTF annotation:

  featureCounts -a hg38_annotation.gtf -o counts.txt subread_results.bam

Tutorials
-------------------
A short tutorial for Subread - http://bioinf.wehi.edu.au/subread
A short tutorial for Subjunc - http://bioinf.wehi.edu.au/subjunc
A short tutorial for featureCounts - http://bioinf.wehi.edu.au/featureCounts
A short tutorial for exactSNP - http://bioinf.wehi.edu.au/exactSNP

Users Guide
--------------
Users Guide can be found in the 'doc' subdirectory of this software package or via URL (http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf).

Citation
--------------
If you use Subread or Subjunc aligners, please cite:

Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013

If you use the featureCounts program, please cite:

Liao Y, Smyth GK and Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014

Mailing lists
--------------
Please post your questions/suggestions to Bioconductor support site(https://support.bioconductor.org/) or Subread google group (https://groups.google.com/forum/#!forum/subread).

