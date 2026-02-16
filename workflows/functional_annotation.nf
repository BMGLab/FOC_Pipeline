// 
include { FILTER_CORE_GENES } from '../modules/filter_core_genes.nf'
include { extractProteins } from '../modules/extract_proteins.nf'

nextflow.enable.dsl = 2

workflow {
    // Read AUGUSTUS and Liftoff outputs directly from the filesystem
    augustus_gff_ch = Channel.fromPath("data/foc/augustus/*/*.gff")
    augustus_faa_ch = Channel.fromPath("data/foc/augustus/*/*.faa")
    liftoff_gff_ch  = Channel.fromPath("data/foc/liftoff/*/*.gff")
    liftoff_fasta_ch = Channel.fromPath("data/foc/liftoff/*/*.fasta")

    // Collect all outputs
    merged_inputs = augustus_gff_ch.combine(augustus_faa_ch)
                                   .combine(liftoff_gff_ch)
                                   .combine(liftoff_fasta_ch)

    // Filter shared core genes (global step)
    core_genes = FILTER_CORE_GENES(merged_inputs)

    core_genes.out_core_filtered
        .filter { it.name.endsWith(".faa") }
        .map { file -> tuple(file.simpleName.replace("_core", ""), file) }
        .set { filtered_faa_by_sample }

    core_genes.out_core_filtered
        .filter { it.name.endsWith(".gff") }
        .map { file -> tuple(file.simpleName.replace("_core", ""), file) }
        .set { filtered_gff_by_sample }


    // AUGUSTUS core genes: antismash, blastp, dbcan
    antismash(augustus_filtered_faa_by_sample)
    blastp(augustus_filtered_faa_by_sample)
    dbcan(augustus_filtered_faa_by_sample)

    // Liftoff core genes: extract proteins → targetp + signalp + wolfpsort
    extracted_proteins = extractProteins(filtered_gff_by_sample)

    targetp(extracted_proteins)
    signalp(extracted_proteins)
    wolfpsort(extracted_proteins)
}