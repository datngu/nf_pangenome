#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-pangenome
========================================================================================
                Pangenome Variant Calling Pipeline with Nextflow.
                [https://github.com/datngu/nf-pangenome](https://github.com/datngu/nf-pangenome)
                Author: Dat T Nguyen
                Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
    Pipeline: Short variant calling with pangenome reference
    
    Workflow:
    1. INDEX_GRAPH       - Build Giraffe indexes for HPRC pangenome graph
    2. GIRAFFE_ALIGN     - Align reads with vg Giraffe (state-of-the-art)
    3. PROJECT_BAM       - Project GAM to BAM format (hg38 coordinates)
    4. CALL_SNPS_INDELS  - Call variants with Pangenome-Aware DeepVariant
    
    Input:
    - HPRC pangenome graph (GBZ format)
    - Reference genome (FASTA, hg38/GRCh38)
    - Sample metadata CSV (sample,read1,read2)
    
    Output:
    - VCF files with SNPs and indels per sample (hg38 coordinates)
    - Optional: BAM alignments (hg38 coordinates)
----------------------------------------------------------------------------------------
*/



/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/test_data/GRCh38.fa"
params.hprc_graph      = "$baseDir/test_data/hprc-v1.1-mc-grch38.gbz"
params.meta_csv        = "$baseDir/samples.csv"
params.outdir          = "test_results"
params.read_type       = "short"  // "short" or "long"
params.output_bam      = false



nextflow.enable.dsl=2



/*
 * INDEX_GRAPH: Build distance and minimizer indexes for Giraffe alignment
 * FIXED: Removed vg autoindex --gbz-input (invalid flag)
 */
process INDEX_GRAPH {
    tag "build_giraffe_indexes"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '128 GB'
    cpus 16

    publishDir "${params.outdir}/indexes", mode: 'copy'
    
    input:
    path graph

    output:
    path "index.gbz", emit: gbz
    path "index.dist", emit: dist
    path "index.min", emit: min
    path "index.zipcodes", emit: zipcodes
    tuple path("index.gbz"), path("index.dist"), path("index.min"), path("index.zipcodes"), emit: indexes
    
    script:
    """
    # Copy input GBZ
    cp ${graph} index.gbz
    
    # Build distance index
    vg index -t ${task.cpus} -j index.dist index.gbz
    
    # Build minimizer index with zipcodes
    vg minimizer -t ${task.cpus} -d index.dist -o index.min -z index.zipcodes index.gbz
    """
}



/*
 * GIRAFFE_ALIGN: Align reads to pangenome graph with vg Giraffe
 */
process GIRAFFE_ALIGN {
    tag "${sample}"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '64 GB'
    cpus 16

    publishDir "${params.outdir}/alignments/${sample}", mode: 'copy', pattern: "*.gam"
    
    input:
    tuple val(sample), path(read1), path(read2)
    tuple path(gbz), path(dist), path(min), path(zipcodes)

    output:
    tuple val(sample), path("${sample}.gam"), emit: gam
    
    script:
    """
    vg giraffe \
        -Z ${gbz} \
        -d ${dist} \
        -m ${min} \
        -z ${zipcodes} \
        -f ${read1} -f ${read2} \
        -t ${task.cpus} \
        --sample ${sample} \
        > ${sample}.gam
    """
}




/*
 * PROJECT_BAM: Project GAM alignment to linear reference and create BAM
 * Note: vg surject projects to GRCh38/hg38 reference paths embedded in the GBZ
 */
process PROJECT_BAM {
    tag "${sample}"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '64 GB'
    cpus 16


    publishDir "${params.outdir}/alignments/${sample}", mode: 'copy', enabled: params.output_bam
    
    input:
    tuple val(sample), path(gam)
    path gbz
    path reference


    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    
    script:
    """
    # Surject to BAM, outputting hg38 coordinates
    vg surject \
        -x ${gbz} \
        -t ${task.cpus} \
        --prune-low-cplx \
        --interleaved \
        --sample ${sample} \
        --read-group "ID:${sample} SM:${sample} LB:lib1 PL:ILLUMINA PU:unit1" \
        --bam-output \
        ${gam} | \
    samtools sort -@ ${task.cpus} -m 4G -o ${sample}.bam -
    
    samtools index -@ ${task.cpus} ${sample}.bam
    """
}



/*
 * CALL_SNPS_INDELS: Call SNPs and Indels with Pangenome-Aware DeepVariant
 */
process CALL_SNPS_INDELS {
    tag "${sample}"
    container 'docker://google/deepvariant:pangenome_aware_deepvariant-1.8.0'
    memory '64 GB'
    cpus 16


    publishDir "${params.outdir}/variants/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(bam), path(bai)
    path reference
    path gbz


    output:
    tuple val(sample), path("${sample}.deepvariant.vcf.gz"), path("${sample}.deepvariant.vcf.gz.tbi"), emit: vcf
    tuple val(sample), path("${sample}.deepvariant.g.vcf.gz"), path("${sample}.deepvariant.g.vcf.gz.tbi"), emit: gvcf
    
    script:
    def model = params.read_type == 'long' ? 'PACBIO' : 'WGS'
    """
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
        --model_type=${model} \
        --ref=${reference} \
        --reads=${bam} \
        --pangenome=${gbz} \
        --output_vcf=${sample}.deepvariant.vcf.gz \
        --output_gvcf=${sample}.deepvariant.g.vcf.gz \
        --num_shards=${task.cpus}
    
    tabix -p vcf ${sample}.deepvariant.vcf.gz
    tabix -p vcf ${sample}.deepvariant.g.vcf.gz
    """
}




workflow {
    // Read samples from CSV
    // CSV format: sample,read1,read2
    samples_ch = Channel
        .fromPath(params.meta_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }


    // Load reference files
    graph_ch = Channel.fromPath(params.hprc_graph)
    genome_ch = Channel.fromPath(params.genome)


    // Index HPRC graph for Giraffe alignment
    INDEX_GRAPH(graph_ch)


    // Align reads to pangenome graph with vg Giraffe
    GIRAFFE_ALIGN(samples_ch, INDEX_GRAPH.out.indexes)


    // Project GAM to BAM (hg38 coordinates)
    PROJECT_BAM(GIRAFFE_ALIGN.out.gam, INDEX_GRAPH.out.gbz, genome_ch)


    // Call variants with Pangenome-Aware DeepVariant
    CALL_SNPS_INDELS(PROJECT_BAM.out.bam, genome_ch, INDEX_GRAPH.out.gbz)
}
