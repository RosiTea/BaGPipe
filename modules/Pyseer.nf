process Pyseer_Unitig {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path unitig
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes ${params.phenotypes} \
        --kmers unitigs.txt.gz \
        --similarity phylogeny_K.tsv \
        --phenotype-column ${params.antibiotic} \
        --output-patterns kmer_patterns_${params.antibiotic}.txt \
        --cpu ${task.cpus} > gwas_${params.antibiotic}_kmers.txt
    count_patterns.py kmer_patterns_${params.antibiotic}.txt > kmer_pattern_count_${params.antibiotic}.txt
    qq_plot.py gwas_${params.antibiotic}_kmers.txt
    # cat <(head -1 gwas_${params.antibiotic}_kmers.txt) <(awk "\$4<1.67E-02 {print \$0}" gwas_${params.antibiotic}_kmers.txt) > significant_variants.txt

    """
}

process Pyseer_Variants {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path variants
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes ${params.phenotypes} \
        --vcf ${variants} \
        --similarity ${k_matrix} \
        --phenotype-column ${params.antibiotic} \
        --output-patterns var_patterns_${params.antibiotic}.txt \
        --cpu ${task.cpus} > gwas_${params.antibiotic}_var.txt
    count_patterns.py var_patterns_${params.antibiotic}.txt > var_pattern_count_${params.antibiotic}.txt
    qq_plot.py gwas_${params.antibiotic}_var.txt
    # cat <(head -1 gwas_${params.antibiotic}_var.txt) <(awk "\$4<1.67E-02 {print \$0}" gwas_${params.antibiotic}_var.txt) > significant_variants.txt

    """
}

process Pyseer_PreAbs {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path pre_abs
    path pheno
    path k_matrix

    output:
    path "*", emit: pyseer_out

    script:
    """
    pyseer \
        --lmm \
        --phenotypes ${params.phenotypes} \
        --pres gene_presence_absence.Rtab \
        --similarity phylogeny_K.tsv \
        --phenotype-column ${params.antibiotic} \
        --output-patterns gene_patterns_${params.antibiotic}.txt \
        --cpu ${task.cpus} > gwas_${params.antibiotic}_preabs.txt
    count_patterns.py gene_patterns_${params.antibiotic}.txt > gene_pattern_count_${params.antibiotic}.txt
    qq_plot.py gwas_${params.antibiotic}_preabs.txt
    # cat <(head -1 gwas.txt) <(awk "\$4<1.67E-02 {print \$0}" gwas.txt) > significant_variants.txt

    # Visualisation
    # cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") \\
    # <(paste <(sed '1d' gwas.txt | cut -d "_" -f 2) \\
    # <(sed '1d' gwas.txt | cut -f 4) | \\
    # awk '{p = -log(\$2)/log(10); print "26",".",\$1,p,p,"0"}' ) | \\
    # tr ' ' '\t' > gwas.plot

    """
}