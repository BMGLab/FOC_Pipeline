// 
process makeBLASTDB {
    publishDir "${params.output_dir}/blast_db", mode: 'copy'
    
    input:
    path reference
    
    output:
    path "blastdb*", emit: blast_db
    
    script:
    """
    source $params.conda_shell
    conda activate blast
    makeblastdb \
        -in ${reference} \
        -dbtype nucl \
        -out blastdb \
        -parse_seqids
    conda deactivate
    """
}