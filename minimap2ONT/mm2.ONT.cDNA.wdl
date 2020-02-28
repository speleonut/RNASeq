workflow minimap2_ONT_cDNA {
    call mimimap2
}
task mimimap2 {
    String htslib
    String samtools
    File program
    File genomeBuild
    File seqFile
    String readGroupID
    String sampleName
    String platform
    String LB
    String outputDir
    Int cores
    command {
        module load ${htslib}
        module load ${samtools}
        ${program} -ax splice \
        -R "@RG\\tID:${readGroupID}\\tLB:${LB}\\tPL:${platform}\\tSM:${sampleName}" \
        -t ${cores} \
        ${refSeq} ${seqFile} |\
        samtools view -bT ${genomeBuild} - |\
        samtools sort -l 5 -m 4G -@${cores} -T${sampleName} -o ${outputDir}/${sampleName}.sort.bam -
        samtools index ${outputDir}/${sampleName}.sort.bam
    }
    output {
        File sortedBAM = "${outputDir}/${sampleName}.sort.bam"
        File sortedBAM = "${outputDir}/${sampleName}.sort.bam.bai"
    }
    runtime {
        job_title: "mm2ont-cDNA"
        slurm_out: "/fast/users/%u/log/mm2ont-cDNA-slurm-%j.out"
        nodes: "1"
        cores: "9"
        queue: "batch"
        requested_memory_mb_per_core: "4"
        time: "05:00:00"
    }
}