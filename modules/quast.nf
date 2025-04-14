//  Quast module for assembly quality assessment
process QUAST {
    tag "$id"
    publishDir "${params.output_dir}/assembly_qc/${id}", mode: 'copy'

    input:
    val id
    path contigs

    output:
    path "quast_output", emit: results

    script:
    """
    source $params.conda_shell
    export JAVA_HOME=\$HOME/miniconda3/envs/quast
    export JAVA_LD_LIBRARY_PATH=\${JAVA_LD_LIBRARY_PATH:-}
    conda activate quast
    quast.py ${contigs} \
             -o quast_output \
             --fast
    conda deactivate
    """
}