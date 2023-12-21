#!/usr/bin/env nextflow

//params.genomes = "$projectDir/genomes/*.fasta"
log.info """\
    FASTA - TO - GWAS PIPELINE
    ==========================
    genomes: ${params.genomes}
    genus: ${params.genus}
    phenotypes: ${params.phenotypes}
    """
    .stripIndent()


// Validation for user-set parameters
//

// Annotate fasta genome assemblies using Prokka
process ProkkaAnnotate {
    publishDir "${params.outdir}/prokka", mode:'copy', overwrite: true

    input:
    path genomes

    output:
    path "${prefix}", emit: prokka_output
    path "${prefix}/*.gff", emit: prokka_output_gff

    script:
    prefix = genomes.getBaseName()
    """
    prokka --cpus ${task.cpus} --genus ${params.genus} --usegenus --outdir ${prefix} --prefix ${prefix} ${genomes}
    """
}

// Build a pangenome, including a multiple sequence alignment of core genes (MAFFT), using Panaroo 
process PanarooAnalysis {
    publishDir "${params.outdir}/panaroo", mode: 'copy', overwrite: true

    input:
    path gff_files

    output:
    path "panaroo_output", emit: panaroo_out
    path "panaroo_output/core_gene_alignment.aln", emit: panaroo_output_core_aln
    path "panaroo_output/gene_presence_absence.Rtab", emit: panaroo_output_pre_abs

    script:
    """
    panaroo -i *.gff -o panaroo_output --clean-mode strict --remove-invalid-genes -a core
    """
}

// Build a phylogeny from the core gene alignment using IQ-TREE
process PhylogeneticAnalysis {
    publishDir "${params.outdir}/iqtree", mode: 'copy', overwrite: true

    input:
    path alignment

    output:
    path "phylogenetic_tree", emit: phylo_tree

    script:
    """
    iqtree -s core_gene_alignment.aln -pre core_tree -nt AUTO -fast -m GTR
    """
}

// Pyseer
process Pyseer {
    publishDir "${params.outdir}/pyseer", mode: 'copy', overwrite: true
    container "quay.io/rositea/tea"

    input:
    path tree 
    path pre_abs_rtab

    output:
    path "pyseer_output", emit: pyseer_output

    script:
    """
    phylogeny_distance.py --lmm ./iqtree_output/core_tree.treefile > pyseer_output/phylogeny_K.tsv
    pyseer --lmm --phenotypes ${params.phenotype} --pres ${pre_abs_rtab} --similarity ./pyseer_output/phylogeny_K.tsv --phenotype-column Tetracycline --output-patterns ./pyseer_output/gene_patterns.txt > ./pyseer_output/gwas.txt
    """
}

workflow {
    genomes_ch = Channel.fromPath(params.genomes)
    genomes_ch.view { "Genomes: $it" }
    ProkkaAnnotate(genomes_ch)
    gff_files = ProkkaAnnotate.out.prokka_output_gff.collect()
    //            .map { dir -> "${dir}/*/*.gff"}
    gff_files.view()
    PanarooAnalysis(gff_files)

    alignment = PanarooAnalysis.out.panaroo_output_core_aln()
    alignment.view()
    PhylogeneticAnalysis(alignment)
}

