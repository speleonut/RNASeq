#!/bin/bash
#SBATCH -J st-final
#SBATCH -o /hpcfs/users/%u/log/stringtie-final-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch            	                                # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8              	                                # number of cores (here uses 8)
#SBATCH --time=06:00:00    	                                # time allocation, which has the format (D-HH:MM), here set to 30 min
#SBATCH --mem=24G         	                                # memory pool for all cores (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                                    # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au                      # Email to which notification will be sent

#Set some paths
HISAT2_INDEXES=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/HISAT2
stringtiePath=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/stringtie_latest
userDir=/hpcfs/users/${USER}
DefaultOutDir=${userDir}/BAM/RNASeq # You might change this depending on how you set up your environment

# StringTie.final.Phoenix.sh
usage()
{
echo "# StringTie.final.Phoenix.sh slurm submission script. Generates read counts for RNASeq analysis
# Before running count number of samples using: wc -l input file
#
# Dependencies:  An input text file 
#
# Usage: sbatch --array 0-(nSamples-1) $0 -f /path/to/inputFile.txt -m /path/to/merged.gtf -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sample ID, you can use your HISAT2 input file listed in the form \"sampleID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -m	OPTIONAL (with caveats). Merged GTF from the previous stringtie --merge.  Will look in \$outDir for the file with defailt name mergedOutput.gtf
# -o	OPTIONAL. Path to where you want to find your files. Default is $DefaultOutDir
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 10/07/2017
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
# 31/08/2020; Mark Corbett; Update for new file paths on Phoenix
#
" 
}

# Parse script options
while [ "$1" != "" ]; do
	case $1 in
        -f )    shift
                seqFile=$1
			    ;;
        -m )    shift
                mergedGTF=$1
                ;;
        -o )    shift
                outDir=$1
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
	echo "#ERROR: You need to specify the path and name of the sample file list"
	exit 1
fi

if [ -z "$outDir" ]; then # If output directory not specified then use the default.
    outDir=$DefaultOutDir
    echo "Using $outDir as the working directory"
fi

if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

if [ -z "$mergedGTF" ]; then # If merged gtf not specified then try to find it
    if [ -s "$outDir/mergedOutput.gtf" ]; then
    	mergedGTF=$outDir/mergedOutput.gtf
    else
    	usage
    	echo "#ERROR: You need to specify -m merged gtf file. I looked for it in $outDir but I couldn't find it"
    	exit 1
    fi
fi
# Define files for the array
sampleID=($(awk -F" " '{print $1}' $seqFile))

# Run final round of StringTie
$stringtiePath/stringtie $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam \
-e -B -G $mergedGTF -p 8 \
-o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.forBallgown.gtf \
-A $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run2.gene_abund.txt \
-C $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run2.cov_refs.gtf
