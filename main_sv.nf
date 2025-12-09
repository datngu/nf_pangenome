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
    1. EXTRACT_REFERENCE    - Extract GRCh38 reference from HPRC graph
    2. DOWNLOAD_HPRC_VCF    - Download pre-processed HPRC VCF for PanGenie
    3. FIX_VCF_CHROMS       - Convert VCF chromosome names to standard notation
    4. PANGENIE_INDEX       - Build PanGenie index from VCF and reference
    5. PANGENIE_GENOTYPE    - Genotype SVs with PanGenie
    
    Input:
    - HPRC pangenome graph (GBZ format) - for reference extraction
    - Sample metadata CSV (sample,reads) - reads can be FASTA/FASTQ
    
    Output:
    - Genotyped VCF files per sample with SV calls
----------------------------------------------------------------------------------------
*/

/*
 Define the default parameters
*/ 
params.hprc_graph      = "$baseDir/test_data/hprc-v1.1-mc-grch38.gbz"
params.hprc_vcf_url    = "https://zenodo.org/records/6797328/files/cactus_filtered_ids.vcf.gz"
params.meta_csv        = "$baseDir/samples.csv"
params.outdir          = "test_results"
params.ref_prefix      = "GRCh38#0#"

nextflow.enable.dsl=2

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
    
    script:
    """"

    vg paths -x ${gbz} -F -S ${params.ref_prefix} > GRCh38.hprc.fa

    sed "s/>${params.ref_prefix}/>/" GRCh38.hprc.fa > GRCh38.fa

    samtools faidx GRCh38.fa

    cut -f1 GRCh38.fa.fai | head -25
    """
}


/*
 * DOWNLOAD_HPRC_VCF: Download pre-processed HPRC VCF for PanGenie
 * NOTE: This VCF is already filtered and formatted for PanGenie use
 */
process DOWNLOAD_HPRC_VCF {
    tag "download_hprc_vcf"
    container 'docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    memory '16 GB'
    cpus 4

    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    output:
    path "hprc.vcf.gz", emit: vcf
    path "hprc.vcf.gz.tbi", emit: tbi
    
    script:
    """
    # Download pre-processed HPRC VCF from Zenodo
    wget -O hprc.vcf.gz ${params.hprc_vcf_url}
    
    # Index VCF
    tabix -p vcf hprc.vcf.gz
    
    # Check VCF contents
    echo "VCF chromosomes:"
    bcftools query -f '%CHROM\\n' hprc.vcf.gz | sort -u | head -10
    
    echo "Total variants:"
    bcftools view -H hprc.vcf.gz | wc -l
    """
}

/*
 * FIX_VCF_CHROMS: Convert VCF chromosome names if needed
 * HPRC VCF may have GRCh38#0# prefix that needs to be removed
 */
process FIX_VCF_CHROMS {
    tag "fix_vcf_chroms"
    container 'docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    memory '32 GB'
    cpus 8

    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    input:
    path vcf
    path tbi
    path reference
    path ref_fai

    output:
    path "hprc.fixed.vcf.gz", emit: vcf
    path "hprc.fixed.vcf.gz.tbi", emit: tbi
    
    script:
    """
    # Check if VCF has GRCh38#0# prefix
    FIRST_CHROM=\$(bcftools query -f '%CHROM\\n' ${vcf} | head -1)
    
    if [[ "\$FIRST_CHROM" == GRCh38#0#* ]]; then
        echo "VCF has HPRC notation (GRCh38#0#), converting to standard (chr)..."
        
        # Create chromosome renaming map
        bcftools view -h ${vcf} | grep "^##contig" | \\
            sed 's/##contig=<ID=//; s/,.*//' | \\
            grep "GRCh38#0#" | \\
            awk '{old=\$1; new=\$1; sub(/^GRCh38#0#/,"",new); print old "\\t" new}' > chr_rename.txt
        
        # Rename chromosomes
        bcftools annotate --rename-chrs chr_rename.txt ${vcf} -Oz -o hprc.fixed.vcf.gz
    else
        echo "VCF already has standard notation, copying..."
        cp ${vcf} hprc.fixed.vcf.gz
    fi
    
    # Index output VCF
    tabix -p vcf hprc.fixed.vcf.gz
    
    # Normalize VCF against reference
    bcftools norm -f ${reference} -Oz -o hprc.fixed.norm.vcf.gz hprc.fixed.vcf.gz
    mv hprc.fixed.norm.vcf.gz hprc.fixed.vcf.gz
    tabix -p vcf hprc.fixed.vcf.gz
    
    echo "Final VCF chromosomes:"
    bcftools query -f '%CHROM\\n' hprc.fixed.vcf.gz | sort -u | head -10
    """
}

/*
 * PANGENIE_INDEX: Build PanGenie index (preprocessing step)
 * This only needs to be done once for all samples
 */
process PANGENIE_INDEX {
    tag "pangenie_index"
    container 'docker://quay.io/biocontainers/pangenie:4.2.1--h7ff8a90_0'
    memory '128 GB'
    cpus 16

    publishDir "${params.outdir}/pangenie_index", mode: 'copy'
    
    input:
    path vcf
    path tbi
    path reference
    path ref_fai

    output:
    path "pangenie-index-*", emit: index_files
    path "pangenie-index.done", emit: done
    
    script:
    """
    # Decompress VCF for PanGenie (it requires uncompressed VCF)
    bcftools view ${vcf} > hprc.vcf
    
    # Run PanGenie-index
    PanGenie-index \
        -v hprc.vcf \
        -r ${reference} \
        -o pangenie-index \
        -t ${task.cpus}
    
    # Create completion marker
    touch pangenie-index.done
    
    # List output files
    echo "PanGenie index files created:"
    ls -lh pangenie-index-*
    """
}

/*
 * PANGENIE_GENOTYPE: Genotype SVs for each sample using PanGenie
 */
process PANGENIE_GENOTYPE {
    tag "${sample}"
    container 'docker://quay.io/biocontainers/pangenie:4.2.1--h7ff8a90_0'
    memory '64 GB'
    cpus 16

    publishDir "${params.outdir}/genotypes/${sample}", mode: 'copy'
    
    input:
    tuple val(sample), path(reads)
    path index_files
    path index_done

    output:
    tuple val(sample), path("${sample}_genotyping.vcf"), emit: vcf
    tuple val(sample), path("${sample}_genotyping.vcf.gz"), path("${sample}_genotyping.vcf.gz.tbi"), emit: vcf_gz
    
    script:
    """
    # Run PanGenie genotyping
    PanGenie \
        -i ${reads} \
        -f pangenie-index \
        -o ${sample} \
        -j ${task.cpus} \
        -t ${task.cpus} \
        -s ${sample}
    
    # Compress and index output VCF
    bgzip -c ${sample}_genotyping.vcf > ${sample}_genotyping.vcf.gz
    tabix -p vcf ${sample}_genotyping.vcf.gz
    
    # Summary statistics
    echo "Genotyping completed for ${sample}"
    echo "Total variants genotyped:"
    grep -v "^#" ${sample}_genotyping.vcf | wc -l
    """
}

workflow {
    // Read samples from CSV
    // CSV format: sample,reads (reads can be FASTA or FASTQ, can be gzipped)
    samples_ch = Channel
        .fromPath(params.meta_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.reads)) }

    samples_ch.view { "Sample: ${it[0]}, Reads: ${it[1]}" }

    // Extract reference from HPRC graph
    EXTRACT_REFERENCE(params.hprc_graph)

    // Download pre-processed HPRC VCF
    DOWNLOAD_HPRC_VCF()

    // Fix VCF chromosome names to match reference
    FIX_VCF_CHROMS(
        DOWNLOAD_HPRC_VCF.out.vcf,
        DOWNLOAD_HPRC_VCF.out.tbi,
        EXTRACT_REFERENCE.out.fasta,
        EXTRACT_REFERENCE.out.fai
    )

    // Build PanGenie index (once for all samples)
    PANGENIE_INDEX(
        FIX_VCF_CHROMS.out.vcf,
        FIX_VCF_CHROMS.out.tbi,
        EXTRACT_REFERENCE.out.fasta,
        EXTRACT_REFERENCE.out.fai
    )

    // Genotype each sample
    PANGENIE_GENOTYPE(
        samples_ch,
        PANGENIE_INDEX.out.index_files.collect(),
        PANGENIE_INDEX.out.done.collect()
    )
}
