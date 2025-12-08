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
    path "hprc_graph.vcf.gz", emit: vcf
    path "hprc_graph.vcf.gz.tbi", emit: vcf_idx
    
    script:
    """
    vg deconstruct -P "GRCh38#0#" -e -a -t ${task.cpus} ${graph} | bgzip > hprc_graph.vcf.gz
    tabix -p vcf hprc_graph.vcf.gz
    """
}

workflow {
    graph_ch = Channel.fromPath(params.hprc_graph)
    EXTRACT_VCF_FROM_GRAPH(graph_ch)
}
