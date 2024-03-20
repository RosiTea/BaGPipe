#!/usr/bin/env nextflow

log.info """\
    FASTA - TO - GWAS PIPELINE
    ==========================
    genus: ${params.genus}
    phenotypes: ${params.phenotypes}
    antibiotic: ${params.antibiotic} 
    reference: ${params.reference} 
    reference_gff: ${params.reference_gff} 
    """
    .stripIndent()

/*
========================================================================================
    HELP
========================================================================================
*/

def printHelp() {
    log.info """
    Usage:
    nextflow run main.nf

    Options:

      --manifest                   Manifest containing paths to FASTQ files (mandatory)
      --phenotypes                 A tab file containing phenotypes for all samples (mandatory)
      --genus                      Genus name for samples (mandatory if starting from FASTA files)
      --genotype_method		   Genotype method to run GWAS, from a choice of three (unitig|pa|snp) (mandatory)
				   Note: unitig is recommended.
      --reference                  Reference genome fasta file (mandatory for significant kmer analysis) 
      --reference_gff              Reference genome gff file (mandatory for significant kmer analysis)
      --mygff                      Input already annotated GFF files; must match sample_ids in manifest (optional)
      --mytree                     Input user preferred phylogenetic tree (optional)
      --pa                         Run GWAS using gene presence and absence, instead of default: unitigs (optional)
      --snp                        Run GWAS using gene variances, instead of default: unitigs (optional)
      --mvcf                       Input already mergerd vcf.gz file (optional)
      --fe                         Run GWAS using fixed model (SEER) (optional)
      --help                       print this help message (optional)

    Alternative Options for Some Processes:

    [PanarooAnalysis]
      --panaroo_clean_mode	   Default: "strict"
      --panaroo_alignment	   Default: "core"
      --panaroo_aligner		   Default: "mafft"
      --panaroo_core_threshold     Default: 0.95
    [PhylogeneticAnalysis]
      --iqtree_model		   Default: "GTR"
    [Pyseer]
      --min_af			   Default: 0.01
      --max_af			   Default: 0.99

    """.stripIndent()
}

if (params.help) {
    printHelp()
    exit(0)
}


/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { validate_parameters } from './modules/helper_functions'
include { ProkkaAnnotate } from './modules/ProkkaAnnotate'
include { PanarooAnalysis } from './modules/PanarooAnalysis'
include { PhylogeneticAnalysis } from './modules/PhylogeneticAnalysis'
include { PyseerKinshipMatrix } from './modules/PyseerKinshipMatrix'
include { UnitigCaller } from './modules/UnitigCaller'
include { MergeVCF } from './modules/MergeVCF'
include { PyseerGenotypeMatrix } from './modules/PyseerGenotypeMatrix'
include { PyseerUnitig; PyseerPreAbs; PyseerVariants } from './modules/Pyseer'
include { SignificantKmers; KmerMap; WriteReferenceText; AnnotateKmers; GeneHitPlot} from './modules/SignificantKmerAnalysis'


/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

validate_parameters()

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    genomes_ch = manifest_ch.splitCsv(header: true, sep: ',')
        .map{ row -> tuple(row.sample_id, row.assembly_path) }

    pheno = Channel.fromPath(params.phenotypes)


    if (params.mytree) {
        tree = Channel.fromPath(params.mytree)

    }
    else {
        if (params.mygff) {
            PanarooAnalysis(gff_files)
            alignment = PanarooAnalysis.out.panaroo_output_core_aln

            PhylogeneticAnalysis(alignment)
            tree = PhylogeneticAnalysis.out.phylo_tree
        }
        else {
            ProkkaAnnotate(genomes_ch)
            gff_files = ProkkaAnnotate.out.prokka_output_gff
                .map { it -> it[1]}
                .collect()
            
            PanarooAnalysis(gff_files)
            alignment = PanarooAnalysis.out.panaroo_output_core_aln
            
            PhylogeneticAnalysis(alignment)
            tree = PhylogeneticAnalysis.out.phylo_tree
        }
    }


    if (genotype_method == "pa"){
        PyseerKinshipMatrix(tree)
        k_matrix = PyseerKinshipMatrix.out.kinship_matrix
        pre_abs = PanarooAnalysis.out.panaroo_output_pre_abs

        PyseerPreAbs(pre_abs,pheno,k_matrix)
    }
    else if (genotype_method == "snp"){
        if (params.mvcf){
            // But currently this has a problem: Must use distance matrix with fixed effects (SEER)
            // This is the classifical way of using MDS 
            // Distance matrix needs mash (1); or we could extract patristic distances from a phylogeny (2). 
                // 1>> mash sketch -s 10000 -o samples *.fa \ mash dist samples.msh samples.msh | square_mash > mash.tsv
                // 2>> python scripts/phylogeny_distance.py core_genome.tree > phylogeny_distances.tsv
            // But, we can use core_gene VCF to generate kinship matrix and use FaST-LMM model?
                // similarity_pyseer --vcf core_gene_snps.vcf sample_list.txt > genotype_kinship.tsv
            if (params.fe){
                // fixed effect way of doing GWAS
            }
            else {
                variants = Channel.fromPath(params.mvcf)
                PyseerGenotypeMatrix(variants, manifest_ch)
                k_matrix = PyseerGenotypeMatrix.out.kinship_matrix
                PyseerVariants(variants,pheno,k_matrix)
            }
            
        }
        else {
            // Either use snippy to call variant from fasta, given a reference, then use Process: MergeVCF
            // Or ask user to input another manifest containing directory of all vcf, then use Process: MergeVCF
        }
    }    
    else if (genotype_method == "unitig"){
        PyseerKinshipMatrix(tree)
        k_matrix = PyseerKinshipMatrix.out.kinship_matrix

        UnitigCaller(manifest_ch)
        unitig = UnitigCaller.out.unitig_out

        PyseerUnitig(unitig,pheno,k_matrix)
        pyseer_result = Pyseer_Unitig.out.pyseer_out

        SignificantKmers(pyseer_result)
        sig_kmer = SignificantKmers.out.sig_kmer_out

        if (params.reference && params.reference_gff){
            ref = Channel.fromPath(params.reference)
            ref_gff = Channel.fromPath(params.reference_gff)

            KmerMap(sig_kmer, ref)
            WriteReferenceText(manifest_ch,ref,gff_files)
            reftxt = WriteReferenceText.out.write_ref_text_out

            AnnotateKmers(sig_kmer,ref,reftxt,ref_gff,gff_files)
            genehit = AnnotateKmers.out.annotated_kmers_out

            GeneHitPlot(genehit)
	}
    }
    else {
	println "Please use a correct genotype method (unitig|pa|snp)."
    }

    
    
}

