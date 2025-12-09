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
    1. EXTRACT_REFERENCE    - Extract GRCh38 reference from HPRC graph (with HPRC notation)
    2. INDEX_GRAPH          - Build Giraffe indexes for HPRC pangenome graph
    3. GIRAFFE_ALIGN        - Align reads with vg Giraffe
    4. PROJECT_BAM          - Project GAM to BAM format (HPRC notation)
    5. CALL_SNPS_INDELS     - Call variants with Pangenome-Aware DeepVariant
    6. FIX_VCF_CHROMS       - Convert chromosome names from HPRC to standard notation
    
    Input:
    - HPRC pangenome graph (GBZ format)
    - Sample metadata CSV (sample,read1,read2)
    
    Output:
    - VCF files with standard chromosome notation (chr1, chr2, etc.)
    - Optional: BAM alignments with HPRC notation
----------------------------------------------------------------------------------------
*/

/*
 Define the default parameters
*/ 
params.hprc_graph      = "$baseDir/test_data/hprc-v1.1-mc-grch38.gbz"
params.meta_csv        = "$baseDir/samples.csv"
params.outdir          = "test_results"
params.read_type       = "short"  // "short" or "long"
params.output_bam      = false
params.ref_prefix      = "GRCh38#0#"  // HPRC v1.1 path prefix

nextflow.enable.dsl=2

/*
 * EXTRACT_REFERENCE: Extract GRCh38 reference from pangenome graph
 * Keeps HPRC notation (GRCh38#0#chr1) for consistency with graph
 */
process EXTRACT_REFERENCE {
    tag "extract_grch38"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '32 GB'
    cpus 8

    publishDir "${params.outdir}/reference", mode: 'copy'
    
    input:
    path gbz

    output:
    path "GRCh38.fa", emit: fasta
    path "GRCh38.fa.fai", emit: fai
    path "chrom_rename.txt", emit: rename_map
    
    script:
    """
    # Extract FASTA with HPRC path names (GRCh38#0#chr1, GRCh38#0#chr2, etc.)
    vg paths -x ${gbz} -F -S GRCh38 > GRCh38.fa
    
    # Index the FASTA
    samtools faidx GRCh38.fa
    
    # Create chromosome renaming map for VCF conversion later
    # Format: old_name new_name (one per line)
    cut -f1 GRCh38.fa.fai | awk -F'#' '{
        if (\$3 ~ /^chr[0-9XY]+\$/) {
            print \$0 "\\t" \$3
        }
    }' > chrom_rename.txt
    
    echo "Chromosome renaming map:"
    head -5 chrom_rename.txt
    
    echo "Extracted chromosomes:"
    cut -f1 GRCh38.fa.fai | head -10
    """
}

/*
 * INDEX_GRAPH: Build distance and minimizer indexes for Giraffe alignment
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
 * PROJECT_BAM: Project GAM to BAM with HPRC chromosome notation
 * BAM will have chromosomes named GRCh38#0#chr1, GRCh38#0#chr2, etc.
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
    path ref_fai

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    
    script:
    """

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
 * CALL_SNPS_INDELS: Call variants with Pangenome-Aware DeepVariant
 * VCF will have HPRC chromosome notation (GRCh38#0#chr1)
 */
process CALL_SNPS_INDELS {
    tag "${sample}"
    container 'docker://google/deepvariant:pangenome_aware_deepvariant-1.8.0'
    memory '64 GB'
    cpus 16

    publishDir "${params.outdir}/variants/${sample}", mode: 'copy', pattern: "*hprc.{vcf.gz,vcf.gz.tbi}"
    
    input:
    tuple val(sample), path(bam), path(bai)
    path reference
    path ref_fai
    path gbz

    output:
    tuple val(sample), path("${sample}.hprc.vcf.gz"), path("${sample}.hprc.vcf.gz.tbi"), emit: vcf
    tuple val(sample), path("${sample}.hprc.g.vcf.gz"), path("${sample}.hprc.g.vcf.gz.tbi"), emit: gvcf
    
    script:
    def model = params.read_type == 'long' ? 'PACBIO' : 'WGS'
    """
    /opt/deepvariant/bin/run_pangenome_aware_deepvariant \
        --model_type=${model} \
        --ref=${reference} \
        --reads=${bam} \
        --pangenome=${gbz} \
        --output_vcf=${sample}.hprc.vcf.gz \
        --output_gvcf=${sample}.hprc.g.vcf.gz \
        --num_shards=${task.cpus}
    
    tabix -p vcf ${sample}.hprc.vcf.gz
    tabix -p vcf ${sample}.hprc.g.vcf.gz
    
    """
}

/*
 * FIX_VCF_CHROMS: Convert chromosome names from HPRC to standard notation
 * GRCh38#0#chr1 -> chr1, GRCh38#0#chr2 -> chr2, etc.
 */
process FIX_VCF_CHROMS {
    tag "${sample}"
    container 'docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    memory '16 GB'
    cpus 4

    publishDir "${params.outdir}/variants/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(vcf), path(tbi)
    tuple val(sample), path(gvcf), path(gtbi)
    path rename_map

    output:
    tuple val(sample), path("${sample}.deepvariant.vcf.gz"), path("${sample}.deepvariant.vcf.gz.tbi"), emit: vcf
    tuple val(sample), path("${sample}.deepvariant.g.vcf.gz"), path("${sample}.deepvariant.g.vcf.gz.tbi"), emit: gvcf
    
    script:
    """
    bcftools annotate \
        --rename-chrs ${rename_map} \
        -Oz -o ${sample}.deepvariant.vcf.gz \
        ${vcf}
    
    tabix -p vcf ${sample}.deepvariant.vcf.gz
    
    bcftools annotate \
        --rename-chrs ${rename_map} \
        -Oz -o ${sample}.deepvariant.g.vcf.gz \
        ${gvcf}
    
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

    samples_ch.view { "Sample channel: ${it[0]}" }

    // Extract reference from graph (HPRC notation)
    EXTRACT_REFERENCE(params.hprc_graph)

    // Index HPRC graph for Giraffe alignment
    INDEX_GRAPH(params.hprc_graph)

    // Align reads to pangenome graph
    GIRAFFE_ALIGN(samples_ch, INDEX_GRAPH.out.indexes.collect())

    // Project GAM to BAM (HPRC notation)
    PROJECT_BAM(
        GIRAFFE_ALIGN.out.gam, 
        INDEX_GRAPH.out.gbz.collect(),
        EXTRACT_REFERENCE.out.fasta.collect(),
        EXTRACT_REFERENCE.out.fai.collect()
    )

    // Call variants with Pangenome-Aware DeepVariant (HPRC notation)
    CALL_SNPS_INDELS(
        PROJECT_BAM.out.bam,
        EXTRACT_REFERENCE.out.fasta.collect(),
        EXTRACT_REFERENCE.out.fai.collect(),
        INDEX_GRAPH.out.gbz.collect()
    )

    // Fix VCF chromosome names to standard notation
    FIX_VCF_CHROMS(
        CALL_SNPS_INDELS.out.vcf,
        CALL_SNPS_INDELS.out.gvcf,
        EXTRACT_REFERENCE.out.rename_map.collect()
    )
}
