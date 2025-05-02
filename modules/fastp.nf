// 
process fastp {
    tag "$sample_id"
    publishDir "${params.output_dir}/fastp/${dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2), val(dir)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html
    path "${sample_id}_fastp.json", emit: json

    script:
    """
    source $params.conda_shell
    conda activate fastp
    fastp -i ${read1} -I ${read2} \\
          -o ${sample_id}_trimmed_R1.fastq.gz \\
          -O ${sample_id}_trimmed_R2.fastq.gz \\
          -h ${sample_id}_fastp.html \\
          -j ${sample_id}_fastp.json
    conda deactivate
    """
}