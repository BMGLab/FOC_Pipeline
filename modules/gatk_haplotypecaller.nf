//
process GATK_Haplotypecaller {
    tag "$sample_id"
    publishDir "${params.output_dir}/gatk", mode: 'copy'

    input:
    val(sample_id)
    path(sorted_bam)
    path(reference)

    output:
    path "${sample_id}.vcf.gz", emit: vcf

    script:
    """
    source $params.conda_shell
    conda activate gatk
    gatk CreateSequenceDictionary \
    -R ../ref/GCF_000149955.1_ASM14995v2_genomic.fna \
    -O ../ref/GCF_000149955.1_ASM14995v2_genomic.dict
    gatk HaplotypeCaller \
        -R $reference \
        -I $sorted_bam \
        -O ${sample_id}.vcf.gz \
        -ERC GVCF
    conda deactivate
    """
}