#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-pangenome
========================================================================================
                Pangenome Variant Calling Pipeline with nextflow.
                [https://github.com/datngu/nf-pangenome](https://github.com/datngu/nf-pangenome)
                Author: Dat T Nguyen
                Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/

/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/test_data/GRCh38.fa"
params.hprc_graph      = "$baseDir/test_data/hprc-v1.1-mc-grch38.gbz"
params.phased_vcf      = "$baseDir/data/1000GP_ONT_shapeit5-phased-callset_final-vcf.phased.vcf.gz"
params.meta_csv        = "$baseDir/samples.csv"
params.outdir          = "test_results"

nextflow.enable.dsl=2

process EXTRACT_VCF_FROM_GRAPH {
    tag "HPRC_graph"
    container 'docker://quay.io/vgteam/vg:v1.65.0'
    memory '128 GB'
    cpus 16

    publishDir "${params.outdir}/graph_vcf", mode: 'copy'
    
    input:
    path graph

    output:
    path "hprc_graph.vcf", emit: vcf
    
    script:
    """
    vg deconstruct -P "GRCh38#0#" -e -a -t ${task.cpus} ${graph} > hprc_graph.vcf
    """
}


process NORMALIZE_VCF {
    tag "normalize_vcf"
    container 'quay.io/biocontainers/bcftools:1.20--h8b25389_0'
    memory '16 GB'
    cpus 4

    publishDir "${params.outdir}/graph_vcf", mode: 'copy'
    
    input:
    path vcf
    path genome

    output:
    path "hprc_graph.normalized.vcf", emit: normalized_vcf
    
    script:
    """
    # Create chromosome renaming map
    bcftools view -h ${vcf} | grep "^##contig" | \\
        sed 's/##contig=<ID=//; s/,.*//' | \\
        grep "GRCh38#0#" | \\
        awk '{print \$1 "\\t" \$1}' | \\
        sed 's/GRCh38#0#//' > chr_rename.txt
    
    # Rename chromosomes in both header and variants
    bcftools annotate --rename-chrs chr_rename.txt ${vcf} -Ov -o renamed.vcf

    # Normalize the VCF  
    bcftools norm -f ${genome} -Ov -o hprc_graph.normalized.vcf renamed.vcf
    """
}



process PANGENIE_INDEX {
    tag "PanGenie_index"
    container 'docker://mgibio/pangenie:v4.2.1-bookworm'
    memory '128 GB'
    cpus 16

    publishDir "${params.outdir}/pangenie_index", mode: 'copy'
    
    input:
    path vcf
    path reference

    output:
    path "pangenie_index*", emit: pangenie_index
    
    script:
    """
    PanGenie-index -v ${vcf} -r ${reference} -o pangenie_index -t ${task.cpus}
    """
}



workflow {
    graph_ch = Channel.fromPath(params.hprc_graph)
    genome_ch = Channel.fromPath(params.genome)

    EXTRACT_VCF_FROM_GRAPH(graph_ch)
    NORMALIZE_VCF(EXTRACT_VCF_FROM_GRAPH.out.vcf, genome_ch)
    PANGENIE_INDEX(NORMALIZE_VCF.out.normalized_vcf, genome_ch)

}
