#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-pangenome-sv
========================================================================================
                Pangenome SV Genotyping Pipeline with PanGenie
                [https://github.com/datngu/nf-pangenome](https://github.com/datngu/nf-pangenome)
                Author: Dat T Nguyen
                Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
    Pipeline: SV genotyping with PanGenie using HPRC pangenome
    
    Workflow:
    1. EXTRACT_REFERENCE      - Extract GENOME_REF reference from HPRC graph
    2. PANGENIE_INDEX         - Build PanGenie index from VCF and reference
    3. PANGENIE_GENOTYPE      - Genotype SVs with PanGenie
    4. CONVERT_TO_BIALLELIC   - Convert multiallelic genotypes to biallelic format
    
    Input:
    - HPRC pangenome graph (GBZ format) - for reference extraction
    - HPRC multiallelic VCF - for variant catalog
    - HPRC biallelic callset VCF - for REF/ALT allele definitions (post-processing)
    - Sample metadata CSV (sample,read1,read2) - reads can be FASTA/FASTQ
    
    Output:
    - Genotyped VCF files per sample with SV calls
----------------------------------------------------------------------------------------
*/


/*
 Define the default parameters
*/ 
params.hprc_graph                = "$baseDir/test_data/hprc-v1.1-mc-grch38.gbz"
params.hprc_pangenie_vcf         = "$baseDir/test_data/pangenie_hprc.vcf"
params.hprc_pangenie_callset_vcf = "$baseDir/test_data/pangenie_hprc_callset.vcf"


params.meta_csv                 = "$baseDir/samples.csv"
params.outdir                   = "test_results"
params.ref_prefix               = "GRCh38"


nextflow.enable.dsl=2


/*
 * EXTRACT_REFERENCE: Extract {GENOME_REF} reference from pangenome graph
 * REMOVE HPRC notation for standardization (chr1, chr2, etc.)
 */

process EXTRACT_REFERENCE {
    tag "extract_${params.ref_prefix}"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '32 GB'
    cpus 8

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path gbz

    output:
    path "genome_ref.fa", emit: fasta
    path "genome_ref.fa.fai", emit: fai
    path "chrom_rename.txt", emit: rename_map

    script:
    """
    # Extract FASTA from the graph
    vg paths -x ${gbz} -F -S ${params.ref_prefix} > genome_ref.fa

    # Index the FASTA
    samtools faidx genome_ref.fa

    # Create chromosome renaming map: old_name -> standard_name
    cut -f1 genome_ref.fa.fai | awk -F'#' '{if (\$3 ~ /^chr[0-9XY]+\$/) print \$0 "\\t" \$3}' > chrom_rename.txt

    # Get HPRC prefix from first line
    PREFIX=\$(cut -f1 chrom_rename.txt | head -1 | sed 's/chr[0-9XY]*\$//')

    # Remove prefix from FASTA and FAI
    sed -i "s/^>\${PREFIX}/>/" genome_ref.fa
    sed -i "s/^\${PREFIX}//" genome_ref.fa.fai
    """
}



/*
 * PANGENIE_INDEX: Build PanGenie index (preprocessing step)
 * This only needs to be done once for all samples
 */
process PANGENIE_INDEX {
    tag "pangenie_index"
    container 'docker://mgibio/pangenie:v4.2.1-bookworm'
    memory '128 GB'
    cpus 16


    publishDir "${params.outdir}/pangenie_index", mode: 'copy'
    
    input:
    path vcf
    path reference
    path ref_fai


    output:
    path "pangenie_index*", emit: index_files
    
    script:
    """
    
    PanGenie-index \
        -v ${vcf} \
        -r ${reference} \
        -o pangenie_index \
        -t ${task.cpus}
    

    """
}


/*
 * PANGENIE_GENOTYPE: Genotype SVs for each sample using PanGenie
 */
process PANGENIE_GENOTYPE {
    tag "${sample}"
    container 'docker://mgibio/pangenie:v4.2.1-bookworm'
    memory '64 GB'
    cpus 16


    publishDir "${params.outdir}/genotypes/", mode: 'copy'
    
    input:
    tuple val(sample), path(read1), path(read2)
    path index_files


    output:


    tuple val(sample), path("${sample}_genotyping.vcf.gz"), path("${sample}_genotyping.vcf.gz.tbi"), emit: vcf_gz
    
    script:
    """
    zcat ${read1} ${read2} > reads.fq


    PanGenie \
        -i reads.fq \
        -f pangenie_index \
        -o ${sample} \
        -j 8 \
        -t 8 \
        -s ${sample}
    
    bgzip -c ${sample}_genotyping.vcf > ${sample}_genotyping.vcf.gz
    tabix -p vcf ${sample}_genotyping.vcf.gz
    
    rm reads.fq ${sample}_genotyping.vcf

    """
}




/*
 * CONVERT_TO_BIALLELIC: Convert multiallelic genotypes to biallelic format
 * Post-processing step for PanGenie output
 * Uses the biallelic callset VCF as reference for REF/ALT alleles
 */
process CONVERT_TO_BIALLELIC {
    tag "${sample}"
    container 'docker://ndatth/ubuntu:22.04'
    memory '16 GB'
    cpus 1

    publishDir "${params.outdir}/biallelic_genotypes/", mode: 'copy'
    
    input:
    tuple val(sample), path(genotype_vcf), path(genotype_tbi)
    path biallelic_vcf

    output:
    tuple val(sample), path("${sample}_genotyping_biallelic.vcf.gz"), path("${sample}_genotyping_biallelic.vcf.gz.tbi"), emit: vcf_gz
    
    script:
    """

    # Convert to biallelic format
    zcat ${genotype_vcf} | \
        python ${projectDir}/bin/convert-to-biallelic.py ${biallelic_vcf} | \
        bgzip -c > ${sample}_genotyping_biallelic.vcf.gz
    

    # Index the output VCF
    tabix -p vcf ${sample}_genotyping_biallelic.vcf.gz
    """
}



workflow {
    // Read samples from CSV
    // CSV format: sample,read1,read2
    samples_ch = Channel
        .fromPath(params.meta_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }


    //samples_ch.view { "Sample channel: ${it[0]}" }


    // Extract reference from HPRC graph
    EXTRACT_REFERENCE(params.hprc_graph)

    // Build PanGenie index (once for all samples)
    PANGENIE_INDEX(
        params.hprc_pangenie_vcf,
        EXTRACT_REFERENCE.out.fasta,
        EXTRACT_REFERENCE.out.fai
    )

    // Genotype each sample
    PANGENIE_GENOTYPE(
        samples_ch,
        PANGENIE_INDEX.out.index_files.collect()
    )

    // Convert genotypes to biallelic format (post-processing)
    CONVERT_TO_BIALLELIC(
        PANGENIE_GENOTYPE.out.vcf_gz,
        params.hprc_pangenie_callset_vcf
    )
}
