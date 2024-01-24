// Annotate fasta genome assemblies using Prokka
process ProkkaAnnotate {
    tag "${sample_id}"
    container "quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_5"
    publishDir "${params.outdir}/prokka", mode:'copy', overwrite: true

    input:
    tuple val(sample_id), path(assembly_path)

    output:
    tuple val(sample_id), path("${sample_id}"), emit: prokka_output
    tuple val(sample_id), path("${sample_id}/*.gff"), emit: prokka_output_gff

    script:
    //prefix = assembly_path.getName()
    """
    prokka \
        --cpus ${task.cpus} \
        --genus ${params.genus} \
        --usegenus \
        --outdir ${sample_id} \
        --prefix ${sample_id} ${assembly_path}
    """
}
