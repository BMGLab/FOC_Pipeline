// 
process extractProteins {
    tag "$id"

    input:
    val id
    path assembly
    path gff

    output:
    val id, emit: sample_id
    path "${id}_predicted_proteins.faa", emit: proteins

    script:
    """
    source $params.conda_shell
    conda activate agat

    agat_sp_fix_cds_phases.pl \
        --gff $gff \
        --fasta ${assembly} \
        -o ${id}_predicted_genes_fixed.gff

    agat_sp_extract_sequences.pl \
        --gff ${id}_predicted_genes_fixed.gff \
        --fasta ${assembly} \
        --type CDS \
        --protein \
        -o ${id}_predicted_proteins.faa

    conda deactivate
    """
}