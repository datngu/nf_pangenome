process AUGMENT_GRAPH {
    tag "augment"
    label 'vg'
    publishDir "${params.outdir}/graphs", mode: 'copy'
    
    input:
    path graph
    path kgp_sv_vcf
    path reference
    
    output:
    path "augmented.gbz", emit: augmented_graph
    
    script:
    """
    vg view -g ${graph} > base.gfa
    vg augment base.gfa ${kgp_sv_vcf} -t ${task.cpus} -a direct -A augmented.gaf > augmented.gfa
    vg gbwt -E -G augmented.gfa -o augmented.gbwt
    vg convert -g augmented.gfa -v > augmented.vg
    vg index -x augmented.xg augmented.vg
    vg gbwt -x augmented.xg -g augmented.gbwt --gbz-format -o augmented.gbz
    rm -f base.gfa augmented.gfa augmented.vg augmented.xg augmented.gbwt augmented.gaf
    """
}
