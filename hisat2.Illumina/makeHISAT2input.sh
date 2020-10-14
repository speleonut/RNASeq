#!/bin/bash

# A short script to make up input for HISAT2
# Usage: makeHISAT2Input.sh /path/to/sequence/files
# 
# Tidied up from a previous Tophat input script
# Mark Corbett 14/10/2020

# Check to see if script has the correct input
if [ -z "$1" ]; then
    seqPath=$(pwd)
    echo "You didn't specify the path to your sequence files so I'm taking a punt that you want to use the current directory a.k.a
$seqPath
If that wasn't what you wanted then run this script again like so:
Usage: $0 /path/to/sequence/files"
else
    seqPath=$1
fi

if [ ! -d "$seqPath" ]; then
   echo "Usage: $0 /path/to/sequence/files
## ERROR: $seqPath is not a directory"
    exit 1
fi

# Remember my current location
currentDir=$(pwd)

# Fetch sample ID from the fastq.gz filenames. YMMV
cd $seqPath
ls *.gz | cut -f1 -d"_" | sort | uniq > $currentDir/tmp.30r4kjfdbit9givids8o093jd.SampleNames.txt # Psuedorandom temporary filename

while read i; do
	ls $seqPath | grep ^$i\_ > $i.files.txt
	sed -i "s,^,$seqPath,g" $i.files.txt	
	grep R1 $i.files.txt > R1.tmp.txt
	grep R2 $i.files.txt > R2.tmp.txt
	rm $i.files.txt
	# for this next bit the code comes from
	# http://unix.stackexchange.com/questions/114943/can-sed-replace-new-line-characters
	# it works!
	sed ':a;N;$!ba;s/\n/,/g' R1.tmp.txt >> R1.tmp.list.txt
	sed ':a;N;$!ba;s/\n/,/g' R2.tmp.txt >> R2.tmp.list.txt
	rm R1.tmp.txt R2.tmp.txt
done < $currentDir/tmp.30r4kjfdbit9givids8o093jd.SampleNames.txt

paste $1 R1.tmp.list.txt R2.tmp.list.txt > $currentDir/RNAList.txt
rm R1.tmp.list.txt R2.tmp.list.txt $currentDir/tmp.30r4kjfdbit9givids8o093jd.SampleNames.txt

cd $currentDir
