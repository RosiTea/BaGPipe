#!/usr/bin/env nextflow

// Pipeline input parameters
params.genomes = "$projectDir/genomes/*.fasta"
params.prokka_output = "$projectDir/prokka_output"
params.panaroo_output = "$projectDir/panaroo_output"
params.iqtree_output = "$projectDir/iqtree_output"
params.pyseer_output = "$projectDir/pyseer_output"
log.info """\
    FASTA - TO - GWAS PIPELINE
    ==========================
    genomes: ${params.genomes}
    """
    .stripIndent()


// Validation for user-set parameters
//

// Annotate fasta genome assemblies using Prokka
process ProkkaAnnotate {
    publishDir params.prokka_output, mode:'copy'

    input:
    path genomes

    output:
    path "${prefix}", emit: prokka_output
    path "${prefix}/*.gff", emit: prokka_output_gff

    script:
    prefix = genomes.getBaseName()
    """
    prokka --cpus 16 --genus Vibrio --usegenus --outdir ${prefix} --prefix ${prefix} ${genomes}
    """
}

// Build a pangenome, including a multiple sequence alignment of core genes (MAFFT), using Panaroo 
process PanarooAnalysis {
    publishDir params.panaroo_output, mode: 'copy'

    input:
    path gff_files

    output:
    path "panaroo_output", emit: panaroo_out

    script:
    """
    panaroo -i *.gff -o panaroo_output --clean-mode strict --remove-invalid-genes --a core
    """
}

workflow {
    genomes_ch = Channel.fromPath(params.genomes)
    ProkkaAnnotate(genomes_ch)
    gff_files = ProkkaAnnotate.out.prokka_output_gff.collect()
    //            .map { dir -> "${dir}/*/*.gff"}
    PanarooAnalysis(gff_files)
    gff_files.view()
}

