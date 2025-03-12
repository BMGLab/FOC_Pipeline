process fastqc {
    tag "$sample_id"
    publishDir "${params.output_dir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "${read1_sample_id}_fastqc.html", emit: html_r1
    path "${read1_sample_id}_fastqc.zip", emit: zip_r1
    path "${read2_sample_id}_fastqc.html", emit: html_r2
    path "${read2_sample_id}_fastqc.zip", emit: zip_r2

    script:
    read1_sample_id = read1.getBaseName().replace('.fq', '')
    read2_sample_id = read2.getBaseName().replace('.fq', '')
    """
    source $params.conda_shell
    export JAVA_HOME=/home/sercanozturk/miniconda3/envs/fastqc
    export JAVA_LD_LIBRARY_PATH=\${JAVA_LD_LIBRARY_PATH:-}
    conda activate fastqc
    fastqc ${read1} ${read2} -o .
    conda deactivate
    """
}