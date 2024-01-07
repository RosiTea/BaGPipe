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
    panaroo -i *.gff -o panaroo_output --clean-mode strict --remove-invalid-genes -a core
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
    iqtree -s core_gene_alignment.aln -pre core_tree -nt AUTO -fast -m GTR
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
    pyseer --lmm --phenotypes ${params.phenotypes} --pres gene_presence_absence.Rtab --similarity phylogeny_K.tsv --phenotype-column ${params.antibiotic} --output-patterns gene_patterns.txt > gwas.txt
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
    pre_abs = PanarooAnalysis.out.panaroo_output_pre_abs
    //pre_abs.view()
    pheno = Channel.fromPath(params.phenotypes)
    Pyseer(pre_abs,pheno,k_matrix)
    
}

