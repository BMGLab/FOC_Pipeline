//  Quast module for assembly quality assessment
process QUAST {
    tag "$id"
    publishDir "${params.output_dir}/assembly_qc/quast/${dir}/${id}", mode: 'copy'

    input:
    val(id)
    path(contigs)
    val(dir)

    output:
    path "basic_stats", emit: basic_stats
    path "icarus_viewers", emit: icarus_viewers
    path "icarus.html", emit: icarus
    path "report.*", emit: report
    path "transposed_report.*", emit: transposed_report

    script:
    """
    source $params.conda_shell
    export JAVA_HOME=\$HOME/miniconda3/envs/quast
    export JAVA_LD_LIBRARY_PATH=\${JAVA_LD_LIBRARY_PATH:-}
    conda activate quast
    quast.py ${contigs} \
        -o .
    conda deactivate
    """
}