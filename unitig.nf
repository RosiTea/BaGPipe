#!/usr/bin/env nextflow

//params.genomes = "$projectDir/genomes/*.fasta"
log.info """\
    FASTA - TO - GWAS PIPELINE
    ==========================
    genus: ${params.genus}
    phenotypes: ${params.phenotypes}
    antibiotic: ${params.antibiotic} 
    """
    .stripIndent()

// Create a variety of parameters options for these scripts

// Validation for user-set parameters
//


// Annotate fasta genome assemblies using Prokka
process ProkkaAnnotate {
    tag "${sample_id}"
    // container "quay.io"
    publishDir "${params.outdir}/prokka", mode:'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly_path)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: prokka_output
    tuple val(sample_id), path("${sample_id}/*.gff"), emit: prokka_output_gff

    script:
    //prefix = assembly_path.getName()
    """
    prokka --cpus ${task.cpus} --genus ${params.genus} --usegenus --outdir ${sample_id} --prefix ${sample_id} ${assembly_path}
    """
}

// Build a pangenome, including a multiple sequence alignment of core genes (MAFFT), using Panaroo 
process PanarooAnalysis {
    publishDir "${params.outdir}/panaroo", mode: 'copy', overwrite: true
    // container "quay.io"

    input:
    path(gff_files)

    output:
    path("panaroo_output"), emit: panaroo_out
    path("panaroo_output/core_gene_alignment.aln"), emit: panaroo_output_core_aln
    path("panaroo_output/gene_presence_absence.Rtab"), emit: panaroo_output_pre_abs

    script:
    """
    panaroo \
        -i *.gff \
        -o panaroo_output \
        --clean-mode strict \
        --remove-invalid-genes \
        -a core \
        -t ${task.cpus}
    """
}

// Build a phylogeny from the core gene alignment using IQ-TREE
process PhylogeneticAnalysis {
    publishDir "${params.outdir}/iqtree", mode: 'copy', overwrite: true
    // container "quay.io"

    input:
    path alignment

    output:
    path "*", emit: iqtree_out
    path "*.treefile", emit: phylo_tree

    script:
    """
    iqtree -s core_gene_alignment.aln -pre core_tree -nt ${task.cpus} -fast -m GTR
    """
}

// Pyseer kinship matrix
process PyseerKinshipMatrix {
    publishDir "${params.outdir}/kinship_matrix", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path tree 

    output:
    path "*K.tsv", emit: kinship_matrix

    script:
    """
    phylogeny_distance.py --lmm core_tree.treefile > phylogeny_K.tsv
    """
}

process UnitigCaller {
    publishDir "${params.outdir}/unitigs", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/unitig-counter:1.1.0--h56fc30b_0"

    input:
    path manifest_ch 

    output:
    path "*", emit: unitig_all
    path "output/*.txt.gz", emit: unitig_out

    script:
    """
    cat ${params.manifest} | tr ',' '\t' > assemblies_list.txt
    unitig-counter -strains assemblies_list.txt -output output -nb-cores ${task.cpus} 
    gzip output/unitigs.txt
    """
}

process Pyseer {
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
    pyseer --lmm --phenotypes ${params.phenotypes} --pres gene_presence_absence.Rtab --similarity phylogeny_K.tsv --phenotype-column ${params.antibiotic} --output-patterns gene_patterns_${params.antibiotic}.txt > gwas.txt
    count_patterns.py gene_patterns_${params.antibiotic}.txt > pattern_count_${params.antibiotic}.txt
    # qq_plot.py gwas.txt
    # cat <(head -1 gwas.txt) <(awk "\$4<1.67E-02 {print \$0}" gwas.txt) > significant_variants.txt

    # Visualisation
    # cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") \\
    # <(paste <(sed '1d' gwas.txt | cut -d "_" -f 2) \\
    # <(sed '1d' gwas.txt | cut -f 4) | \\
    # awk '{p = -log(\$2)/log(10); print "26",".",\$1,p,p,"0"}' ) | \\
    # tr ' ' '\t' > gwas.plot

    """
}

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
    pyseer --lmm --phenotypes ${params.phenotypes} --kmers unitigs.txt.gz --similarity phylogeny_K.tsv --phenotype-column ${params.antibiotic} --output-patterns kmer_patterns_${params.antibiotic}.txt --cpu 8 > gwas_${params.antibiotic}_kmers.txt
    count_patterns.py kmer_patterns_${params.antibiotic}.txt > kmer_pattern_count_${params.antibiotic}.txt
    qq_plot.py gwas_${params.antibiotic}_kmers.txt
    # cat <(head -1 gwas_${params.antibiotic}_kmers.txt) <(awk "\$4<1.67E-02 {print \$0}" gwas_${params.antibiotic}_kmers.txt) > significant_variants.txt

    """
}


// 

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    genomes_ch = manifest_ch.splitCsv(header: true, sep: ',')
        .map{ row -> tuple(row.sample_id, row.assembly_path) }

    //genomes_ch.view()
    ProkkaAnnotate(genomes_ch)
    gff_files = ProkkaAnnotate.out.prokka_output_gff
        .map { it -> it[1]}
        .collect()
    //gff_files.view()
    PanarooAnalysis(gff_files)

    alignment = PanarooAnalysis.out.panaroo_output_core_aln
    //alignment.view()
    PhylogeneticAnalysis(alignment)

    tree = PhylogeneticAnalysis.out.phylo_tree
    //tree.view()
    PyseerKinshipMatrix(tree)

    k_matrix = PyseerKinshipMatrix.out.kinship_matrix
    //k_matrix.view()
    //pre_abs = PanarooAnalysis.out.panaroo_output_pre_abs
    //pre_abs.view()
    pheno = Channel.fromPath(params.phenotypes)

    UnitigCaller(manifest_ch)
    unitig = UnitigCaller.out.unitig_out
    unitig.view()

    Pyseer_Unitig(unitig,pheno,k_matrix)
    
}

