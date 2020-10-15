#!/bin/bash
#SBATCH -J hisat2-cDNA
#SBATCH -o /hpcfs/users/%u/log/hisat2-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8              	                                # number of cores (here uses 8)
#SBATCH --time=12:00:00    	                                # time allocation, which has the format (D-HH:MM), here set to 30 min
#SBATCH --mem=24G         	                                # memory pool for all cores (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

#Set some paths
HISAT2_INDEXES=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/HISAT2
hisat2Path=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/hisat2-2.2.1
stringtiePath=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/stringtie_latest
userDir=/hpcfs/users/${USER}

# hisat2.PE.Phoenix.sh
usage()
{
echo "# hisat2.PE.Phoenix.sh slurm submission script for mapping stranded Illumina paired end RNA-seq reads to a genome. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#
# Usage: sbatch --array 0-(n Samples-1) $0 -f inputFile.txt -g GenomeIndex  -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -g	OPTIONAL. Genome index for hisat2 e.g. default is grch38_tran. Check the $HISAT2_INDEXES folder for options.
# -o	OPTIONAL. Path to where you want to find your files default is $userDir/BAM/RNASeq
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
                SeqFile=$1
			    ;;
        -g )    shift
                refSeq=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -h | --help )    usage
                         exit 1
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
    # -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\""
	exit 1
fi


if [ -z "$refSeq" ]; then # If genome not specified then default to grch38_tran
        refSeq=grch38_tran
		echo "# WARN: No genome build selected so defulting to $refSeq"
fi

if [ -z "$outDir" ]; then #If output directory not specified then run in current directory
    outDir=$userDir/BAM/RNASeq
	echo "# INFO: Using $outDir as the working directory"
fi

if [ ! -d "$outDir" ]; then
		mkdir -p $outDir
fi

if [ ! -d "$outDir/logs" ]; then
		mkdir -p $outDir/logs
fi

# Define files for the array
sampleID=($(awk -F" " '{print $1}' $SeqFile))

# Load modules
module load GCC
module load HTSlib
module load SAMtools
module load Python/3.7.0

# Run HISAT2
if [ -d "$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    rm -r $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/*
    echo "INFO: Contents of $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]} cleared!"
else
    mkdir -p $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi
${hisat2Path}/hisat2 -x $HISAT2_INDEXES/$refSeq/$refSeq \
-1 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $2}') \
-2 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $3}') \
--novel-splicesite-outfile $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.novel.ss.txt \
--dta -p 8 \
--new-summary \
| samtools view -b - | samtools sort -@ 8 -o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam
samtools index $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam

# Run first round of StringTie
${stringtiePath}/stringtie $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam \
-G $HISAT2_INDEXES/$refSeq/$refSeq.gff -p 8 \
-o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gtf \
-A $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gene_abund.txt \
-C $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.cov_refs.gtf
