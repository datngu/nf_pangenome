#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.meta_csv = "samples.csv"
params.hprc_graph = "samples.csv"

process INDEX_GRAPH {
    tag "build_giraffe_indexes"
    
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
    # Simulate index creation
    touch index.gbz
    touch index.dist
    touch index.min
    touch index.zipcodes
    
    echo "Created indexes for ${graph}"
    """
}

process GIRAFFE_ALIGN {
    tag "${sample}"
    
    input:
    tuple val(sample), path(read1), path(read2)
    tuple path(gbz), path(dist), path(min), path(zipcodes)
    
    output:
    tuple val(sample), path("${sample}.gam")
    
    script:
    """
    echo "Aligning ${sample}" > ${sample}.gam
    """
}

workflow {
    // Read samples from CSV
    samples_ch = Channel
        .fromPath(params.meta_csv)
        .splitCsv(header: true, strip: true)
        .map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }
    
    samples_ch.view { "Sample channel: ${it[0]}" }
    
    // Load graph
    graph_ch = Channel.fromPath(params.hprc_graph)
    
    // Index graph (runs once)
    INDEX_GRAPH(graph_ch)
    
    INDEX_GRAPH.out.indexes.view { "Index output: ${it}" }
    
    // Align reads (should run twice - once per sample)
    GIRAFFE_ALIGN(samples_ch, INDEX_GRAPH.out.indexes.collect())
    
    // Count final outputs
    GIRAFFE_ALIGN.out.collect().view { gams ->
        println "\n=== GIRAFFE_ALIGN SUMMARY ==="
        println "Total alignments: ${gams.size()}"
    }
}
