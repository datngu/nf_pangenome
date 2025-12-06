process INDEX_GRAPH {
    tag "index"
    label 'vg'
    memory '128 GB'
    cpus 16
    publishDir "${params.outdir}/indexes", mode: 'copy'
    
    input:
    path graph
    
    output:
    path "graph.gbz", emit: gbz
    path "graph.dist", emit: dist
    path "graph.min", emit: min
    tuple path("graph.gbz"), path("graph.dist"), path("graph.min"), emit: indexes
    path "graph.vcf", emit: vcf_uncompressed
    
    script:
    """
    cp ${graph} graph.gbz
    vg index -j graph.dist -t ${task.cpus} graph.gbz
    vg minimizer -k 29 -w 11 -t ${task.cpus} -g graph.gbz -d graph.dist -o graph.min
    vg deconstruct -P GRCh38 -e -a -t ${task.cpus} graph.gbz > graph.vcf
    """
}

process INDEX_GRAPH_VCF {
    tag "compress_vcf"
    label 'bcftools'
    publishDir "${params.outdir}/indexes", mode: 'copy'
    
    input:
    path vcf
    
    output:
    path "graph.vcf.gz", emit: vcf
    path "graph.vcf.gz.tbi", emit: vcf_idx
    
    script:
    """
    bgzip -c ${vcf} > graph.vcf.gz
    tabix -p vcf graph.vcf.gz
    """
}
