#!/bin/bash
#SBATCH -J prepDE
#SBATCH -o /hpcfs/users/%u/log/prepDE-slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 2              	                                # number of cores
#SBATCH --time=01:30:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=8G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# hisat2.PE.Phoenix.sh
#Set some paths and modules
stringtiePath=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/stringtie_latest
userDir=/hpcfs/users/${USER}
module purge
module use /apps/skl/modules/all
modList=("Python/3.7.0" "HTSlib/1.17-GCC-11.2.0")

usage()
{
echo "# hisat2.PE.Phoenix.sh slurm submission script for mapping stranded or unstranded Illumina paired end RNA-seq reads to a genome. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#
# Usage: sbatch $0 -f inputFile.txt -o /path/to/outDir | [-h | --help]
#
# Options: 
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -o	REQUIRED. Path to where you want to find your files default is $userDir/BAM/RNASeq
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
    echo "# ERROR: You need to specify the path and name of the sequence file list
# -f	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\""
    exit 1
fi
if [ -z "$outDir" ]; then #If output directory not specified then run in current directory
    usage
    echo "# ERROR: You need to specify the path to where your counts data sample folders are.
# -o	REQUIRED. path to where your counts data sample folders are."
    exit 1
fi

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

# Run prepDE
paste <(cut -f1 $seqFile) <(find $outDir/*/*.forBallgown.gtf) > $outDir/gtfDE.list.txt

cd $outDir
$stringtiePath/prepDE.py3 -i $outDir/gtfDE.list.txt -l 150

# Tidy up and compress some files
(grep ^"#" $outDir/mergedOutput.gtf; grep -v ^"#" $outDir/mergedOutput.gtf | sort -k1,1 -k4,4n) | bgzip > $outDir/mergedOutput.gtf.gz
tabix $outDir/mergedOutput.gtf.gz
if [ -f "$outDir/mergedOutput.gtf.gz.tbi" ]; then
    rm $outDir/mergedOutput.gtf
fi

find $outDir/*/*.gtf > $outDir/gtf.compress.list.txt
while read gtf ; do
    (grep ^"#" ${gtf}; grep -v ^"#" ${gtf} | sort -k1,1 -k4,4n) | bgzip > ${gtf}.gz
    tabix ${gtf}.gz
done < $outDir/gtf.compress.list.txt

while read gtf ; do
    if [ -f "${gtf}.gz.tbi" ]; then
       rm ${gtf}
    fi
done < $outDir/gtf.compress.list.txt

rm $outDir/gtf.compress.list.txt
