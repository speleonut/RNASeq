#!/bin/bash
# hisat2.launcher.sh
#Set some paths and modules
userDir=/hpcfs/users/${USER}
scriptDir=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/mark/RNASeq/hisat2.Illumina

usage()
{
echo "# hisat2.launcher.sh script that coordinates slurm submission for mapping and counting stranded or unstranded Illumina paired end RNA-seq reads to a genome. 
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#
# Usage: hisat2.launcher.sh $0 -f inputFile.txt -g GenomeIndex  -S RF|FR -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -g	REQUIRED. Genome index for hisat2 use grch38_tran for human grcm38_tran form mouse. Check the $HISAT2_INDEXES folder for other options.
# -S    OPTIONAL (with caveats). Default not set (unstranded). Use either RF or FR depending on your library. You probably want to use RF.
#                                See https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/ for advice you can BLAT a few reads to work it out.
# -o	OPTIONAL. Path to where you want to find your files default is $userDir/BAM/RNASeq
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 20/01/2022
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
#
" 
}

# Parse script options
while [ "$1" != "" ]; do
	case $1 in
        -f )    shift
                seqFile=$1
			    ;;
        -g )    shift
                refSeq=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -S )    shift
                strandFlag=$1
                ;;
        -h | --help )    usage
                         exit 0
                         ;;
        * )    usage
               exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "$seqFile" ]; then #If sequence file list in a text file is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the sequence file list.
    # -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\""
	exit 1
fi
if [ -z "$refSeq" ]; then # If genome not specified then default to grch38_tran
    usage
    echo "# ERROR: No genome build selected."
    echo "# -g	REQUIRED. Genome index for hisat2 use grch38_tran for human grcm38_tran form mouse. Check the $HISAT2_INDEXES folder for other options."
    exit 1
fi
if [ -z "$outDir" ]; then #If output directory not specified then run in current directory
    outDir=$userDir/BAM/RNASeq
fi
if [ -z "$strandFlag" ]; then
    echo "# INFO: Using unstranded mapping protocol."
else
    strandSet="-S $strandFlag"
fi
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi
echo "# INFO: Using $outDir as the working directory."

# Ensure log directory exists
if [ ! -d "$userDir/log" ]; then
    mkdir -p $userDir/log
    echo "# INFO: Your slurm .out files will be found in $userDir/log"
fi

# Set sample number for array jobs
nSamples=$(sed '1d' $seqFile | wc -l)

# Load modules
mapJob=`sbatch --array 0-$nSamples --export=ALL $scriptDir/hisat2.PE.more_mismatches.Phoenix.sh -f $seqFile -g $refSeq -o $outDir $strandSet`
mapJob=$(echo $mapJob | cut -d" " -f4)
mergeJob=`sbatch --export=ALL --dependency=afterok:$mapJob $scriptDir/StringTieMerge.Phoenix.sh -g $refSeq -o $outDir`
mergeJob=$(echo $mergeJob | cut -d" " -f4)
countsJob=`sbatch --array 0-$nSamples --export=ALL --dependency=afterok:$mergeJob $scriptDir/StringTie.final.Phoenix.sh -f $seqFile -o $outDir`
countsJob=$(echo $countsJob | cut -d" " -f4)
sbatch --export=ALL --dependency=afterok:$countsJob $scriptDir/prepDE.sh -f $seqFile -o $outDir
