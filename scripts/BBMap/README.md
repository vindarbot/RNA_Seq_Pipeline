# BBMap/BBTools

(Not Offical) BBMap short read aligner for DNA/RNAseq, and other bioinformatic tools.  
BBTools bioinformatics tools, including BBMap.

I have moved those dozens of shell scripts from root to `./sh/` to make it tidy.

* [SEQanswers Page](http://seqanswers.com/forums/showthread.php?t=41057)
* [SourceForge Page](https://sourceforge.net/projects/bbmap/)

BBMap/BBTools are now open source. Please try it out - it's a 3MB download, and written in pure Java, so installation is trivial - just unzip and run. Handles all sequencing platforms (Illumina, PacBio, 454, Sanger, Nanopore, etc) except Solid colorspace, which I removed to simplify the code.

A Powerpoint comparison of performance (speed, memory, sensitivity, specificity) on various genomes, compared to bwa, bowtie2, gsnap, smalt:

<https://drive.google.com/file/d/0B3llHR93L14wbks0ZURFcFhFR1E/edit?usp=sharing>

...but in summary, BBMap is similar in speed to bwa, with much better sensitivity and specificity than any other aligner I've compared it to. It uses more memory than Burrows-Wheeler-based aligners, but in exchange, the indexing speed is many times faster.

## Compiling

````bash
brew cask install java
cd jni/
JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home make -f makefile.osx
cd -

brew install --with-java open-mpi
ant -Dmpijar=/usr/local/lib/mpi.jar
````

To build Eclipse project with `ant`, you need to put `ecj-4.5.2.jar` to `/usr/local/Cellar/ant/1.9.6/libexec/lib/`.

## Running

To run with jar:

````bash
java -cp ./dist/lib/BBTools.jar align2.BBMap
java -cp ./dist/lib/BBTools.jar align2.BBMap usejni=t
````

There is documentation in the docs folder and displayed by shellscripts when run with no arguments. But for example:

`bbmap.sh ref=ecoli.fa`  
...will build an index and write it to the present directory

`bbmap.sh in=reads.fq out=mapped.sam`  
...will map to the indexed reference

`bbmap.sh in1=reads1.fq in2=reads2.fq out=mapped.sam ref=ecoli.fa nodisk`  
...will build an index in memory and map paired reads to it in a single command

If your OS does not support shellscripts, replace 'bbmap.sh' like this:  
`java -Xmx23g -cp /path/to/current align2.BBMap in=reads.fq out=mapped.sam`

...where /path/to/current is the location of the 'current' directory, and -Xmx23g specifies the amount of memory to use. This should be set to about 85% of physical memory (the symbols 'm' or 'g' specify megs or gigs), or more, depending on your virtual memory configuration. Human reference requires around 21 GB; generally, references need around 7 bytes per base pair, and a minimum of 500 MB at default settings. However, there is a reduced memory mode ('usemodulo') that only needs half as much memory. The shellscripts are just wrappers that display usage information and set the -Xmx parameter. 

## Original Readme

BBMap/BBTools readme  
Written by Brian Bushnell  
Last updated December 23, 2015  

The BBTools package was written by Brian Bushnell, with the exception of the (optional, but faster) C, JNI, and MPI components, which were written by Jonathan Rood.

All tools in the BBTools package are free to use.  If you use BBTools in work leading to a publication, and BBTools has not yet been published, please cite it something like this:  
BBMap - Bushnell B. - sourceforge.net/projects/bbmap/

License:

The BBMap package is open source and free to use with no restrictions.  For more information, please read Legal.txt and license.txt.

Documentation:

Documentation is in the /bbmap/docs/ directory, and in each tool's shellscript in /bbmap/.  
readme.txt: This file.  
UsageGuide.txt: Contains basic installation and usage information.  Please read this first!  
ToolDescriptions.txt: Contains a list of all BBTools, a description of what they do, and their hardware requirements.  
compiling.txt: Information on compiling JNI code.  
readme_config.txt: Usage information about config files.  
readme_filetypes.txt: More detailed information on file formats supported by BBTools.  
changelog.txt: List of changes by version, and current known issues.

Tool-specific Guides:

Some tools have specific guides, like BBDukGuide.txt.  They are in /bbmap/docs/guides/.  For complete documentation of a tool, I recommend that you read UsageGuide.txt first (which covers the shared functionality of all tools), then the tool's specific guide if it has one (such as ReformatGuide.txt), then the tool's shellscript (such as reformat.sh) which lists all of the flags.

If you have any questions not answered in the documentation, please look at the relevant SeqAnswers thread (linked from here: http://seqanswers.com/forums/showthread.php?t=41057) and post a question there if it is not already answered.  You can also contact JGI's BBTools team at bbtools@lbl.gov, or me at bbushnell@lbl.gov.  But please read the documentation first.

Special thanks for help with shellscripts goes to:  
Alex Copeland (JGI), Douglas Jacobsen (JGI/NERSC), Bill Andreopoulos (JGI), sdriscoll (SeqAnswers), Jon Rood (JGI/NERSC), and Elmar Pruesse (UC Denver).

Special thanks for helping to support BBTools goes to Genomax (SeqAnswers).
