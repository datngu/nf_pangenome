process CALL_SNPS_INDELS {
    tag "${sample}"
    label 'deepvariant'
    publishDir "${params.outdir}/variants/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path reference
    path gbz
    val read_type
    
    output:
    tuple val(sample), path("${sample}.deepvariant.vcf.gz"), path("${sample}.deepvariant.vcf.gz.tbi"), emit: vcf
    tuple val(sample), path("${sample}.deepvariant.g.vcf.gz"), emit: gvcf
    
    script:
    def model = read_type == 'long' ? 'PACBIO' : 'WGS'
    """
    mkdir -p intermediate logs
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
        --model_type=${model} \
        --ref=${reference} \
        --reads=${bam} \
        --pangenome=${gbz} \
        --output_vcf=${sample}.deepvariant.vcf.gz \
        --output_gvcf=${sample}.deepvariant.g.vcf.gz \
        --num_shards=${task.cpus} \
        --intermediate_results_dir=./intermediate \
        --logging_dir=./logs \
        --make_examples_extra_args="min_mapping_quality=1"
    tabix -p vcf ${sample}.deepvariant.vcf.gz
    tabix -p vcf ${sample}.deepvariant.g.vcf.gz
    rm -rf intermediate
    """
}
