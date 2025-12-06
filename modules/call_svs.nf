process CALL_SVS {
    tag "${sample}"
    label 'pangenie'
    publishDir "${params.outdir}/variants/${sample}", mode: 'copy', pattern: "*.pangenie_genotyping.vcf"
    
    input:
    tuple val(sample), path(read1), path(read2)
    path vcf
    path vcf_idx
    path reference
    
    output:
    tuple val(sample), path("${sample}.pangenie_genotyping.vcf"), emit: vcf
    
    script:
    """
    if [[ ${vcf} == *.gz ]]; then
        gunzip -c ${vcf} > pangenie_input.vcf
    else
        cp ${vcf} pangenie_input.vcf
    fi
    
    cat ${read1} ${read2} > combined_reads.fastq.gz
    
    PanGenie-index -v pangenie_input.vcf -r ${reference} -o pangenie_index -t ${task.cpus}
    PanGenie -f pangenie_index -i combined_reads.fastq.gz -s ${sample} -j ${task.cpus} -t ${task.cpus} -o pangenie_output
    PanGenie-vcf -f pangenie_index -z pangenie_output_genotyping.cereal -s ${sample} -o ${sample}.pangenie
    
    rm -f pangenie_input.vcf combined_reads.fastq.gz
    rm -rf pangenie_index pangenie_output*
    """
}

process COMPRESS_PANGENIE_VCF {
    tag "${sample}"
    label 'bcftools'
    publishDir "${params.outdir}/variants/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf)
    
    output:
    tuple val(sample), path("${sample}.pangenie.vcf.gz"), emit: vcf
    
    script:
    """
    bgzip -c ${vcf} > ${sample}.pangenie.vcf.gz
    tabix -p vcf ${sample}.pangenie.vcf.gz
    """
}
