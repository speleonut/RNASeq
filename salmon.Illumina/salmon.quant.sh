#!/bin/bash
#SBATCH -J fishyQuant
#SBATCH -o /hpcfs/users/%u/log/salmon_quant-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 12              	                                # number of cores
#SBATCH --time=02:00:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=36G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# salmon.quant.sh
#Set some paths and define functions
salmon_prog=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/salmon_latest/bin/salmon
salmon_index_dir=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/salmon
annotation_dir=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/annotations
userDir=/hpcfs/users/${USER}
threads=12 # Set the same as n above

# Genome list (alter case statement below to add new options)
set_genome_build() {
case "${buildID}" in
    hg38 | GRCh38 )    buildID="hg38"
                       genomeBuild="$salmon_index_dir/hg38/default"
                       tx2gene="$annotation_dir/Homo_sapiens.hg38.refgenie.txp2gene.tsv"
                       ;;
    mm10 | GRCm38 )    buildID="mm10"
                       genomeBuild="$salmon_index_dir/mm10/default"
                       tx2gene="$annotation_dir/Mus_musculus.mm10.refgenie.txp2gene.tsv"
                       ;;
    * )         echo "## ERROR: Genome build ${buildID} not recognized. Available options are hg38 or mm10."
                exit 1
                ;;
esac
}

usage()
{
echo "# salmon.quant.sh slurm submission script for quantifying transcripts from stranded Illumina paired end RNA-seq reads. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#
# Usage: sbatch --array 0-(n Samples-1) $0 -f inputFile.txt -g GenomeBuild -l library-type -o /path/to/outDir -x \"other salmon flags\" | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -g	REQUIRED. Genome index for salmon e.g. default is grch38_tran. Use either hg38 or mm10 or add to $salmon_index_dir folder and edit this script for new options.
# -l    OPTIONAL (with caveats). The default is to select ISR which should be correct for most Illumina stranded library preparations. 
# -o	OPTIONAL. Path to where you want to find your files default is $userDir/RNASeq/salmon/genomeBuild
# -x    OPTIONAL. Any other salmon program flags you want to set in quoted text (YMMV).
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 17/01/2022
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
#
" 
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -f )    shift
                SeqFile=$1
                ;;
        -g )    shift
                buildID=$1
                ;;
        -l )    shift
                libType=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -x )    shift
                otherFlags=$1
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
if [ -z "$SeqFile" ]; then #If sequence file list in a text file is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the sequence file list
    # -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"sample-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\""
	exit 1
fi
if [ -z "$buildID" ]; then
    usage
    echo "# ERROR: You need to select either hg38 or mm10 as the genome build."
    echo "# -g	REQUIRED. Genome index for salmon e.g. default is grch38_tran. Use either hg38 or mm10 or add to $salmon_index_dir folder and edit this script for new options."
    exit 1
else
    set_genome_build
    echo "# INFO: using $buildID as the genome build code."
fi
if [ -z "$outDir" ]; then #If output directory not specified then run in current directory
    outDir=$userDir/RNASeq/salmon/$buildID
	echo "# INFO: Using $outDir as the working directory"
fi
if [ -z "$libType" ]; then #If sequence file list in a text file is not supplied then do not proceed
    libType="ISR"
    echo "# WARN: Using the default library type ISR.  This is usually what you want for stranded Illumina RNA library preparations.  If you're not sure, map with HISAT2 and check strandedness using IGV."
fi

# Define files for the array
sampleID=($(awk -F" " '{print $1}' $SeqFile))
read1=($(awk -F" " '{print $2}' $SeqFile))
read2=($(awk -F" " '{print $3}' $SeqFile))


if [ ! -d "$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    mkdir -p $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi

# Run salmon
$salmon_prog quant -l $libType \
-i $genomeBuild \
-g $tx2gene \
-1 $(echo ${read1[$SLURM_ARRAY_TASK_ID]} | tr "," " ") \
-2 $(echo ${read2[$SLURM_ARRAY_TASK_ID]} | tr "," " ") \
-p $threads \
-o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]} \
$otherFlags
