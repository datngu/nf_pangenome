process POSTPROCESS {
    tag "${sample}"
    label 'bcftools'
    publishDir "${params.outdir}/final/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(dv_vcf), path(dv_tbi)
    tuple val(sample), path(pangenie_vcf)
    
    output:
    tuple val(sample), path("${sample}.merged_variants.vcf.gz"), path("${sample}.merged_variants.vcf.gz.tbi"), emit: vcf
    tuple val(sample), path("${sample}.snps_indels.vcf.gz"), emit: snps_indels
    tuple val(sample), path("${sample}.svs.vcf.gz"), emit: svs
    
    script:
    """
    bcftools view -i 'STRLEN(REF)<50 && STRLEN(ALT)<50' ${dv_vcf} -Oz -o ${sample}.snps_indels.vcf.gz
    bcftools index -t ${sample}.snps_indels.vcf.gz
    
    cp ${pangenie_vcf} ${sample}.svs.vcf.gz
    tabix -p vcf ${sample}.svs.vcf.gz
    
    bcftools concat -a -Oz -o ${sample}.merged_unsorted.vcf.gz \
        ${sample}.snps_indels.vcf.gz ${sample}.svs.vcf.gz
    bcftools sort ${sample}.merged_unsorted.vcf.gz -Oz -o ${sample}.merged_variants.vcf.gz
    tabix -p vcf ${sample}.merged_variants.vcf.gz
    
    rm -f ${sample}.merged_unsorted.vcf.gz
    """
}
