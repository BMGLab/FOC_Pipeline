process LAST {
    tag "${id}"
    publishDir "${params.output_dir}/assembly_qc/last/${id}", mode: 'copy'

    input:
    val id
    path query
    path reference

    output:
    path "filtered_reads.*", emit: filtered_reads
    path "long_alignments.tab", emit: long_alignments
    path "high_quality_alignments.tab", emit: high_quality_alignments

    script:
    """
    source ${params.conda_shell}
    conda activate last
    # Step 1: Filter query and reference sequences
    seqkit grep -v -r -p "Chr0|NW" ${query} > filtered_query.fasta
    seqkit grep -v -r -p "NW" ${reference} > filtered.fasta

    # Step 2: Create LAST database
    lastdb -P8 reference_db filtered.fasta

    # Step 3: Perform alignment and post-processing
    lastal -P8 reference_db filtered_query.fasta > query.maf
    last-split query.maf > filtered_reads.maf

    # Step 4: Convert MAF to tab and filter alignments
    maf-convert tab query.maf | awk '\$4 > 500' > long_alignments.tab
    maf-convert tab query.maf | awk '\$10 > 90' > high_quality_alignments.tab

    # Step 5: Generate dotplot
    last-dotplot filtered_reads.maf filtered_reads.png
    conda deactivate
    """
}