#!/bin/bash
#SBATCH -J st-merge
#SBATCH -o /hpcfs/users/%u/log/stringtie-merge-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 8              	                                # number of cores
#SBATCH --time=00:30:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=4G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# CHANGE this if you aren't me! Email to which notification will be sent

#Set some paths
HISAT2_INDEXES=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/HISAT2
stringtiePath=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/stringtie_latest
userDir=/hpcfs/users/${USER}
DefaultOutDir=${userDir}/BAM/RNASeq # You might change this depending on how you set up your environment

# StringTieMerge.Phoenix.sh
usage()
{
echo "# StringTieMerge.Phoenix.sh slurm submission script. Feed it a text file to find GFF files called by the hisat2.PE.Phoenix.sh script
#
# Dependencies:  An input text file 
#
# Usage: sbatch $0 -f inputFile.txt -g GenomeIndex  -o /path/to/outDir -p prefixForMerged.gtf | [-h | --help]
#
# Options: 
# -f	OPTIONAL (with caveats). Path and file name of a text file with sequences listed in the form /path/to/out.1.gtf to final gtf with a newline character
#       If not specified the script will try to find the files it needs but success is not guaranteed.
# -g	OPTIONAL. Genome index for hisat2. The default is grch38_tran
# -p    OPTIONAL. Specific prefix for the merged gtf. Default is mergedOutput.gtf
# -o	OPTIONAL. Path to where you want to find your files default is $DefaultOutDir
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 12/07/2017
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
# 31/08/2020; Mark Corbett; Update for new file paths on Phoenix
# 20/01/2022; Mark Corbett; Update annotation file type
#
" 
}

# Parse script options
while [ "$1" != "" ]; do
	case $1 in
        -f )    shift
                gtfList=$1
			    ;;
        -g )    shift
                refSeq=$1
                ;;
        -p )    shift
                gtfPrefix=$1
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


if [ -z "$refSeq" ]; then # If genome not specified then default to grch38_tran
        refSeq=grch38_tran
		echo "# WARN: No genome build selected so defulting to $refSeq"
fi

if [ -z "$gtfPrefix" ]; then #If gtfPrefx not specified then use a generic one
    gtfPrefix="mergedOutput.gtf"
	echo "INFO: Using $gtfPrefix as the output file name"
fi

if [ -z "$outDir" ]; then #If output directory not specified then run in default directory
    outDir=$DefaultOutDir
	echo "INFO: Using $outDir as the output directory"
fi

if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

if [ -z "$gtfList" ]; then #If a gtf file list text file is not supplied then try to find them or die after this horrible nested if sequence!
	find $outDir/*/*.stringtie.run1.gtf > $outDir/tmp.gtfList.txt
	if [ -s "$outDir/tmp.gtfList.txt" ]; then
		gtfList=$outDir/tmp.gtfList.txt
	else
	    find $DefaultoutDir/*/*.stringtie.run1.gtf > $outDir/tmp.gtfList.txt
	    if [ -s "$outDir/tmp.gtfList.txt" ]; then
	    	gtfList=$outDir/tmp.gtfList.txt
	    else
            usage
            echo "#ERROR: You need to specify the path and name of the GTF file list I looked in subfolders of $outDir and $DefaultoutDir to create the file but the GTF can't be located"
            exit 1
        fi
    fi
fi

# Run the script
$stringtiePath/stringtie --merge -p8 -G $HISAT2_INDEXES/$refSeq/$refSeq.gtf \
-o $outDir/$gtfPrefix \
$gtfList
