#!/bin/bash
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8              	                                # number of cores (here uses 8)
#SBATCH --time=12:00:00    	                                # time allocation, which has the format (D-HH:MM), here set to 30 min
#SBATCH --mem=24G         	                                # memory pool for all cores (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=mark.corbett@adelaide.edu.au  	# CHANGE this if you aren't me! Email to which notification will be sent

#Set some paths
HISAT2_INDEXES=/data/neurogenetics/RefSeq/HISAT2

# hisat2.PE.Phoenix.sh
usage()
{
echo "# hisat2.PE.Phoenix.sh slurm submission script feed it a text file formatted to find fastq files for your RNASeq analysis
# Before running count number of samples using: wc -l input file
#
# Dependencies:  An input text file 
#
# Usage: sbatch --array 0-(nSamples-1) $0 -f inputFile.txt -g GenomeIndex  -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -g	REQUIRED. Genome index for hisat2
# -o	OPTIONAL. Path to where you want to find your files default is $FASTDIR/BAM/RNASeq
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 10/07/2017
# email: mark.corbett@adelaide.edu.au
# Modified (Date, Name, Description):
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
                OutDir=$1
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
	echo "#ERROR: You need to specify the path and name of the sequence file list"
	exit 1
fi


if [ -z "$refSeq" ]; then # If genome not specified then do not proceed
	usage
	echo "#ERROR: You need to specify -g genome build
	# -g	REQUIRED. Genome build your exome is mapped to."
 	exit 1
fi

if [ -z "$OutDir" ]; then #If output directory not specified then run in current directory
    OutDir=$FASTDIR/BAM/RNASeq
	echo "Using $OutDir as the working directory"
fi

if [ ! -d $OutDir ]; then
		mkdir -p $OutDir
fi

if [ ! -d $OutDir/logs ]; then
		mkdir -p $OutDir/logs
fi

# Define files for the array
sampleID=($(awk -F" " '{print $1}' $SeqFile))

# Load modules
module load HISAT2/2.1.0-foss-2016b
module load HTSlib/1.3.1-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load StringTie/1.3.3-foss-2017a

# Run HISAT2
if [ -d $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]} ]; then
    rm -r $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/*
    echo "INFO: Contents of $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]} cleared!"
else
    mkdir -p $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi
hisat2 -x $HISAT2_INDEXES/$refSeq/$refSeq \
-1 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $2}') \
-2 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $3}') \
--novel-splicesite-outfile $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.novel.ss.txt \
--dta --rf -p 8 \
--new-summary \
| samtools view -b - | samtools sort -@ 8 -o $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam
samtools index $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam

# Run first round of StringTie
stringtie $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam \
-G $HISAT2_INDEXES/$refSeq/$refSeq.gff -p 8 --rf \
-o $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gtf \
-A $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gene_abund.txt \
-C $OutDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.cov_refs.gtf \

