process SignificantKmers {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true

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

process WriteReferenceText {
    publishDir "${params.outdir}/significant_unitigs", mode: 'copy', overwrite: true

    input:
    path manifest_ch 
    path ref
    path gff_files

    output:
    path "*", emit: write_ref_text_out

    script:
    """
    reference="${params.reference}"
    reference_base=\${reference%.fa}
    echo -n "" > references.txt
    echo "\${reference_base}.fa    \${reference_base}.gff    ref" >> references.txt
    while IFS=, read -r sample_id assembly_path
    do
        if [[ \$sample_id != "sample_id" ]]; then
            echo "\$assembly_path    \$sample_id.gff    draft" >> references.txt
        fi
    done < ${params.manifest}
    """
}

process AnnotateKmers {
    publishDir "${params.outdir}/annotated_unitigs", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path sig_kmer
    path ref
    path reftxt
    path ref_gff
    path gff_files
    

    output:
    path "*", emit: annotated_kmers_out

    script:
    """
    annotate_hits_pyseer ${sig_kmer} ${reftxt} annotated_kmers.txt
    summarise_annotations.py annotated_kmers.txt > gene_hits.tsv
    """
}

process GeneHitPlot {
    publishDir "${params.outdir}/annotated_unitigs", mode: 'copy', overwrite: true
    // container "quay.io/biocontainers/r-ggplot2:2.2.1--r3.3.2_0"
    // It seems that r-ggrepel container alone is sufficient for this R script. 
    container "quay.io/biocontainers/r-ggrepel:0.6.5--r3.3.2_0"

    input:
    path genehit

    output:
    path "*", emit: genehit_plot_out

    script:
    """
    Rscript ${projectDir}/scripts/gene_hit_summary_plot.R 
    """
}
