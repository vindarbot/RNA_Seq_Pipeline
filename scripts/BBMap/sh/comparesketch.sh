#!/bin/bash
#comparesketch in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified July 7, 2016

Description:  Compares a MinHashSketch others, and prints their kmer identity.
Specifically, the first file will be compared to all others.


Usage:  comparesketch.sh in=<file,file,file...>
Alternative:  comparesketch.sh file file file


Standard parameters:
in=<file,file...>   List of sketches.
size=10000          Default size of sketches.  Does not need to be specified,
                    but setting it correctly may let the program run faster.

Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx200m"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
}
calcXmx "$@"

comparesketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
	fi
	local CMD="java $EA $z -cp $CP sketch.SketchTool $@"
	echo $CMD >&2
	eval $CMD
}

comparesketch "$@"
