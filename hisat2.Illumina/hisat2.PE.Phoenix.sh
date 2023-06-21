#!/bin/bash
#SBATCH -J hisat2-cDNA
#SBATCH -o /hpcfs/users/%u/log/hisat2-slurm-%j.out
#SBATCH -p skylake,icelake,v100cpu            	            # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 8              	                                # number of cores
#SBATCH --time=12:00:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=24G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# hisat2.PE.Phoenix.sh
#Set some paths and modules
HISAT2_INDEXES=/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/HISAT2
hisat2Path=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/hisat2-2.2.1
stringtiePath=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/stringtie_latest
userDir=/hpcfs/users/${USER}
module purge
module use /apps/skl/modules/all
modList=("arch/skylake" "GCC/7.3.0" "HTSlib/1.9" "SAMtools/1.8" "Python/3.7.0")

set_strand()
{
case $strandFlag in
    "RF" | "rf" )    hisat_string="--rna-strandness RF --fr" # Yes that isn't a typo
                     st_string="--rf"
                     ;;
    "FR" | "fr" )    hisat_string="--rna-strandness FR --fr"
                     st_string="--fr"
                     ;;
    * )              echo "# WARN: You used $strandFlag to set the strandedness option but this script only accepts RF or FR (not case sensitive).
Mapping will proceed without any strandedness options. You'll still get reads mapped and counts for transcrpts but it might not be what you wanted."
                     ;;
esac
}

usage()
{
echo "# hisat2.PE.Phoenix.sh slurm submission script for mapping stranded or unstranded Illumina paired end RNA-seq reads to a genome. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#
# Usage: sbatch --array 0-(n Samples-1) $0 -f inputFile.txt -g GenomeIndex  -S RF|FR -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -g	OPTIONAL (with caveats). Genome index for hisat2 e.g. default is grch38_tran. Check the $HISAT2_INDEXES folder for options.
# -S    OPTIONAL (with caveats). Default not set (unstranded). Use either RF or FR depending on your library. You probably want to use RF.
#                                See https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/ for advice you can BLAT a few reads to work it out.
# -o	OPTIONAL. Path to where you want to find your files default is $userDir/BAM/RNASeq
# -h | --help	Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 10/07/2017
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
# 31/08/2020; Mark Corbett; Update for new file paths on Phoenix
# 20/01/2022; Mark Corbett; Explicit setting of modules. Add in strandedness flag. Add in strandedness setting. Update annotation to gencode GTF. Deprecate hisat2.PE.stranded.Phoenix.sh
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
fi
if [ ! -z "$strandFlag" ]; then
    set_strand
else
    echo "# INFO: Using unstranded mapping protocol."
fi
if [ ! -d "$outDir" ]; then
		mkdir -p $outDir
fi
echo "# INFO: Using $outDir as the working directory"

if [ ! -d "$outDir/logs" ]; then
		mkdir -p $outDir/logs
fi

# Define list of sample ID for the array and create output
sampleID=($(awk -F" " '{print $1}' $SeqFile))
if [ -d "$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    if [ -f "$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam" ]; then
        mv $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]} $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}.old
        mkdir -p $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
        mv $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}.old $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/
        echo "#WARN: You are re-running the HISAT pipeline, the contents of $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]} have been moved to $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.old to avoid overwrite errors."
    fi
else
    mkdir -p $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

# Run HISAT2
${hisat2Path}/hisat2 -x $HISAT2_INDEXES/$refSeq/$refSeq \
-1 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $2}') \
-2 $(grep ${sampleID[$SLURM_ARRAY_TASK_ID]} $SeqFile | awk -F" " '{print $3}') \
--novel-splicesite-outfile $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.novel.ss.txt \
--dta $hisat_string -p 8 \
--new-summary \
| samtools view -b - | samtools sort -@ 8 -o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam
samtools index $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam

# Run first round of StringTie
${stringtiePath}/stringtie $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.hisat2.bam \
-G $HISAT2_INDEXES/$refSeq/$refSeq.gtf $st_string -p 8 \
-o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gtf \
-A $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.gene_abund.txt \
-C $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.stringtie.run1.cov_refs.gtf
