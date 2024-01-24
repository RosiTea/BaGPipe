process MergeVCF {
    publishDir "${params.outdir}/vcf", mode: 'copy', overwrite: true
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_0"

    input:
    path vcfs 

    output:
    path "*", emit: merged_vcf

    script:
    """
    bcftools merge -m none -0 -O z *.vcf.gz > merged.vcf.gz
    #  Only one ALT variant per row is supported...
    """
}