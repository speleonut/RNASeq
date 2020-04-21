#!/bin/bash
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8              	                                # number of cores (here uses 8)
#SBATCH --time=00:30:00    	                                # time allocation, which has the format (D-HH:MM), here set to 30 min
#SBATCH --mem=4G         	                                # memory pool for all cores (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=mark.corbett@adelaide.edu.au  	# CHANGE this if you aren't me! Email to which notification will be sent

#Set some paths
HISAT2_INDEXES=/data/neurogenetics/RefSeq/HISAT2
DefaultOutDir=$FASTDIR/BAM/RNASeq # You might change this depending on how you set up your \$FASTDIR

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
#       If not specified the script will try to find the files it needs; success is not guaranteed.
# -g	REQUIRED. Genome index for hisat2. eg. grch38_tran
# -p    OPTIONAL. Specific prefix for the merged gtf. Default is mergedOutput.gtf
# -o	OPTIONAL. Path to where you want to find your files default is $FASTDIR/BAM/RNASeq
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 12/07/2017
# email: mark.corbett@adelaide.edu.au
# Modified (Date, Name, Description):
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
                OutDir=$1
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


if [ -z "$refSeq" ]; then # If genome not specified then do not proceed
	usage
	echo "#ERROR: You need to specify -g genome build
	# -g	REQUIRED. Genome build your RNASeq is mapped to."
 	exit 1
fi

if [ -z "$gtfPrefix" ]; then #If gtfPrefx not specified then use a generic one
    gtfPrefix="mergedOutput.gtf"
	echo "INFO: Using $gtfPrefix as the output file name"
fi

if [ -z "$OutDir" ]; then #If output directory not specified then run in default directory
    OutDir=$DefaultOutDir
	echo "INFO: Using $OutDir as the output directory"
fi

if [ ! -d "$OutDir" ]; then
		mkdir -p $OutDir
fi

if [ -z "$gtfList" ]; then #If a gtf file list text file is not supplied then try to find them or die after this horrible nested if sequence!
	find $OutDir/*/*.stringtie.run1.gtf > $OutDir/tmp.gtfList.txt
	if [ -s $OutDir/tmp.gtfList.txt ]; then
		gtfList=$OutDir/tmp.gtfList.txt
	else
	    find $DefaultOutDir/*/*.stringtie.run1.gtf > $OutDir/tmp.gtfList.txt
	    if [ -s "$OutDir/tmp.gtfList.txt" ]; then
	    	gtfList=$OutDir/tmp.gtfList.txt
	    else
            usage
            echo "#ERROR: You need to specify the path and name of the GTF file list I looked in subfolders of $OutDir and $DefaultOutDir to create the file but the GTF can't be located"
            exit 1
        fi
    fi
fi

# Load the modules
module load StringTie/1.3.3-foss-2017a

# Run the script
stringtie --merge -p8 -G $HISAT2_INDEXES/$refSeq/$refSeq.gff \
-o $OutDir/$gtfPrefix \
$gtfList
