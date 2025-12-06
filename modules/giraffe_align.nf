process GIRAFFE_ALIGN {
    tag "${sample}"
    label 'vg'
    publishDir "${params.outdir}/alignments/${sample}", mode: 'copy', enabled: params.output_bam
    
    input:
    tuple val(sample), path(read1), path(read2)
    tuple path(gbz), path(dist), path(min)
    
    output:
    tuple val(sample), path("${sample}.gam"), emit: gam
    tuple val(sample), path(read1), path(read2), emit: reads
    
    script:
    """
    vg giraffe -Z ${gbz} -d ${dist} -m ${min} \
        -f ${read1} -f ${read2} \
        -t ${task.cpus} --sample ${sample} > ${sample}.gam
    """
}
