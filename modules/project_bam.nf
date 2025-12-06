process PROJECT_BAM {
    tag "${sample}"
    label 'vg'
    publishDir "${params.outdir}/alignments/${sample}", mode: 'copy', enabled: params.output_bam
    
    input:
    tuple val(sample), path(gam)
    path graph
    path reference
    
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    tuple val(sample), path(gam), emit: reads
    
    script:
    """
    vg surject -x ${graph} -b -i -t ${task.cpus} \
        --prune-low-cplx --interleaved --sample ${sample} \
        --read-group "ID:${sample} SM:${sample} PL:ILLUMINA" \
        ${gam} > ${sample}.unsorted.bam
    samtools sort -@ ${task.cpus} -m 4G -o ${sample}.bam ${sample}.unsorted.bam
    samtools index -@ ${task.cpus} ${sample}.bam
    rm -f ${sample}.unsorted.bam
    """
}
