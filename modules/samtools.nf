//
process samTools {
    tag "$id"
    publishDir "${params.output_dir}/samtools/${id}", mode: 'copy'

    input:
    tuple val(id), path(sam)

    output:
    tuple val(id), path("${id}_sorted.bam"), path("${id}_sorted.bam.bai"), emit: bam_bai

    script:
    """
    source $params.conda_shell
    conda activate samtools
    samtools view -bS ${sam} | samtools sort -o ${id}_sorted.bam
    samtools index ${id}_sorted.bam
    conda deactivate
    """
}