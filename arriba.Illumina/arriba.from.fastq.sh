#!/bin/bash

#SBATCH -J HeyAmigo
#SBATCH -o /hpcfs/users/%u/log/Arriba-slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 10              	                                # number of cores
#SBATCH --time=12:00:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=160G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# arriba.from.fastq.sh
#Set some paths and define functions
STAR_prog="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/STAR_2.7.11b/Linux_x86_64_static/STAR"
arriba_prog_dir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/arriba_v2.4.0"
RefDir="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq"
STAR_index_dir="$RefDir/STAR"
annotation_dir="$RefDir/annotations"
userDir="/hpcfs/users/${USER}"
threads=8 # Set one less than n above

# Genome list (alter case statement below to add new options)
set_genome_build() {
case "${buildID}" in
    hg38 | GRCh38 )    buildID="GRCh38"  # Chromosomes are 1-22,X,Y with no chr prefix
                       genomeBuild="$RefDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                       GTF="$annotation_dir/Homo_sapiens.GRCh38.111.gtf.gz"
                       blacklist="$arriba_prog_dir/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
                       known_fusions="$arriba_prog_dir/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
                       GFF="$arriba_prog_dir/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
                       ;;
    mm10 | GRCm38 )    buildID="mm10"
                       genomeBuild="$RefDir/GRCm38_68.fa"
                       GTF="$annotation_dir/Mus_musculus.GRCm38.gencode.vM25.chr_patch_hapl_scaff.annotation.gtf"
                       blacklist="$arriba_prog_dir/database/blacklist_mm10_GRCm38_v2.4.0.tsv.gz"
                       known_fusions="$arriba_prog_dir/database/known_fusions_mm10_GRCm38_v2.4.0.tsv.gz"
                       GFF="$arriba_prog_dir/database/protein_domains_mm10_GRCm38_v2.4.0.gff3"
                       ;;
    * )         echo "## ERROR: Genome build ${buildID} not recognized. Available options are GRCh38 or mm10."
                exit 1
                ;;
esac
}

usage()
{
echo "# arriba.from.fastq.sh a slurm submission script for identifying fusion transcripts from Illumina paired end RNA-seq reads. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\"
#                The SLURM log directory must exist ${userDir}/log or submission to SLURM will fail
#
# Usage: sbatch --array 0-(n Samples-1) $0 -f inputFile.txt [ -g GenomeBuild -o /path/to/outDir ] | $0 [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\"
#                 The optional SV file can be a VCF or a tab-delimited file as specified in the ARRIBA documents: https://arriba.readthedocs.io/en/latest/input-files/#structural-variant-calls-from-wgs
# -g	OPTIONAL. Default hg38. Genome index for STAR and Arriba. Use either GRCh38 or mm10 or add a folder with the build ID to $STAR_index_dir folder and edit this script for new options.
# -o	OPTIONAL. Path to where you want to find your files, the default is $userDir/RNASeq/arriba/genomeBuild/. 
#                 Each analyses will be put in a subfolder of this output directory using the sampleID if you use the default, or specifiy the directory yourself.
# -h | --help     Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 11/04/2024
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
if [ -z "$SeqFile" ]; then #If sequence file list in a text file is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the sequence file list
    # -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"sample-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\""
	exit 1
fi
if [ -z "$buildID" ]; then
    buildID="GRCh38"
    echo "# WARN: The default genome build code $buildID was selected."
    set_genome_build
else
    set_genome_build
    echo "# INFO: Using $buildID as the genome build code."
fi
if [ -z "$outDir" ]; then #If output directory not specified then run in current directory
    outDir="$userDir/RNASeq/arriba/$buildID"
	echo "# INFO: Using $outDir as the working directory"
fi

# Define variables for the array jobs
sampleID=($(awk -F" " '{print $1}' $SeqFile))
read1=($(awk -F" " '{print $2}' $SeqFile))
read2=($(awk -F" " '{print $3}' $SeqFile))
SVfile=($(awk -F" " '{print $4}' $SeqFile))

if [ ! -d "$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    mkdir -p $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi

tmpDir="$outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}_STARtmp"
if [ -d "$tmpDir" ]; then
    rm -r $tmpDir  # STAR likes to start with a clean tmp directory
    echo "# WARN: A pre-existing STAR temporary directory $tmpDir was removed before starting this run.  This is usually occurs when a previous run failed."
fi

if [ ! -z "${SVfile[SLURM_ARRAY_TASKID]}" ]; then
    SVparams="-d ${SVfile[SLURM_ARRAY_TASKID]}"
fi

# Do the thing!
$STAR_prog \
    --runThreadN $threads \
    --genomeDir $STAR_index_dir/$buildID --genomeLoad NoSharedMemory --outTmpDir $tmpDir\
    --readFilesIn ${read1[$SLURM_ARRAY_TASK_ID]} ${read2[$SLURM_ARRAY_TASK_ID]} --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
    --outFileNamePrefix $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/ |
$arriba_prog_dir/arriba \
    -x /dev/stdin \
    -o $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.$buildID.fusions.tsv \
    -O $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.$buildID.fusions.discarded.tsv \
    -a $genomeBuild -g $GTF $SVparams \
    -b $blacklist -k $known_fusions -t $known_fusions -p $GFF
    