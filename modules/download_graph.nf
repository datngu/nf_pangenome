process DOWNLOAD_GRAPH {
    tag "HPRC_graph"
    label 'vg'
    memory '4 GB'
    publishDir "${params.outdir}/graphs", mode: 'copy'
    
    output:
    path "hprc.gbz", emit: graph
    
    script:
    """
    # Using curl as it's more commonly available in containers
    curl -L -o hprc.gbz ${params.hprc_graph_url}
    vg stats -z hprc.gbz
    """
}
