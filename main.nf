#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DOWNLOAD_GRAPH } from './modules/download_graph'
include { AUGMENT_GRAPH } from './modules/augment_graph'
include { INDEX_GRAPH; INDEX_GRAPH_VCF } from './modules/index_graph'
include { GIRAFFE_ALIGN } from './modules/giraffe_align'
include { PROJECT_BAM } from './modules/project_bam'
include { CALL_SNPS_INDELS } from './modules/call_snps_indels'
include { CALL_SVS; COMPRESS_PANGENIE_VCF } from './modules/call_svs'
include { POSTPROCESS } from './modules/postprocess'

workflow {
    // Read samples from CSV file
    // CSV format: sample,read1,read2
    samples_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }
    
    reference_ch = Channel.fromPath(params.reference)
    
    // Get or download HPRC graph
    if (params.hprc_graph) {
        graph_ch = Channel.fromPath(params.hprc_graph)
    } else {
        DOWNLOAD_GRAPH()
        graph_ch = DOWNLOAD_GRAPH.out.graph
    }
    
    // Optional: augment graph with 1KGP SVs
    if (params.augment_graph && params.kgp_sv_vcf) {
        kgp_sv_ch = Channel.fromPath(params.kgp_sv_vcf)
        AUGMENT_GRAPH(graph_ch, kgp_sv_ch, reference_ch)
        final_graph_ch = AUGMENT_GRAPH.out.augmented_graph
    } else {
        final_graph_ch = graph_ch
    }
    
    // Index graph (once for all samples)
    INDEX_GRAPH(final_graph_ch)
    INDEX_GRAPH_VCF(INDEX_GRAPH.out.vcf_uncompressed)
    
    // Process each sample
    GIRAFFE_ALIGN(samples_ch, INDEX_GRAPH.out.indexes)
    PROJECT_BAM(GIRAFFE_ALIGN.out.gam, final_graph_ch, reference_ch)
    CALL_SNPS_INDELS(PROJECT_BAM.out.bam, reference_ch, INDEX_GRAPH.out.gbz, params.read_type)
    CALL_SVS(PROJECT_BAM.out.reads, INDEX_GRAPH_VCF.out.vcf, INDEX_GRAPH_VCF.out.vcf_idx, reference_ch)
    COMPRESS_PANGENIE_VCF(CALL_SVS.out.vcf)
    POSTPROCESS(CALL_SNPS_INDELS.out.vcf, COMPRESS_PANGENIE_VCF.out.vcf)
}
