process SignificantKmers {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true
    //container "quay.io/rositea/tea"

    input:
    path pyseer_result

    output:
    path "*", emit: sig_kmer_out

    script:
    """
    threshold=\$(grep 'Threshold:' kmer_pattern_count_${params.antibiotic}.txt | cut -f2)
    cat <(head -1 gwas_${params.antibiotic}_kmers.txt) \
        <(awk -v thresh=\$threshold '\$4<thresh {print \$0}' \
        gwas_${params.antibiotic}_kmers.txt) > significant_kmers.txt
    """
}

process KmerMap {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path sig_kmer
    path ref

    output:
    path "*.plot", emit: kmer_map_out

    script:
    """
    phandango_mapper significant_kmers.txt ${params.reference} mapped_kmers.plot
    """
}