include { fastp } from '../modules/fastp.nf'
include { samTools } from '../modules/samtools.nf'

process Bowtie2 {
    tag "$id"
    publishDir "${params.output_dir}/alignment/${id}", mode: 'copy'

    input:
    tuple val(id), path(read1), path(read2)
    path reference

    output:
    tuple val(id), path("*.sam"), emit: sam

    script:
    """
    source $params.conda_shell
    conda activate bowtie2
    bowtie2-build ${reference} ${id}_bowtie2_index
    bowtie2 -x ${id}_bowtie2_index \
            -1 ${read1} -2 ${read2} \
            -S ${id}_aligned.sam \
            2> ${id}_bowtie2.log
    conda deactivate
    """
}

process BCFtools {
    tag "$id"
    publishDir "${params.output_dir}/bcftools/${id}", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)
    path reference

    output:
    tuple val(id), path("${id}_filtered_variants.vcf"), emit: vcf

    script:
    """
    source $params.conda_shell
    conda activate bcftools
    bcftools mpileup -f ${reference} --max-depth 1000 ${bam} -Ou | bcftools call -mv -Ov -o ${id}_variants.bcf

    bcftools filter \
        -s LowQual \
        -e 'QUAL<20 || INFO/DP<10' \
        -o ${id}_filtered_variants.bcf \
        ${id}_variants.bcf

    bcftools view ${id}_filtered_variants.bcf -Ov -o ${id}_filtered_variants.vcf
    conda deactivate
    """
}

process SnpEff {
    tag "$id"
    publishDir "${params.output_dir}/snpeff/${id}", mode: 'copy'
    
    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}_annotated.vcf"), emit: annotated_vcf
    path "${id}_snpEff_summary.html", emit: summary_html
    path "snpEff_genes.txt", emit: gene_stats
    path "snpEff_summary.csv", emit: summary_csv

    script:
    """
    source $params.conda_shell
    export JAVA_HOME=${params.home}/miniconda3/envs/snpeff
    export JAVA_LD_LIBRARY_PATH=\${JAVA_LD_LIBRARY_PATH:-}
    conda activate snpeff
    export SNPEFF_JAR=${params.home}/miniconda3/envs/snpeff/share/snpeff-5.2-1/snpEff.jar
    java -jar \$SNPEFF_JAR databases | grep ${params.snpeff_db}

    java -Xmx8g -jar \$SNPEFF_JAR \
        -v ${params.snpeff_db} \
        -stats ${id}_snpEff_summary.html \
        -csvStats snpEff_summary.csv \
        ${vcf} > ${id}_annotated.vcf

    grep -v "^#" ${id}_annotated.vcf | cut -f 8 | grep -o 'ANN=.*' | \
    sed 's/ANN=//g' | tr ',' '\\n' | cut -f 4 -d '|' | sort | uniq -c > snpEff_genes.txt
    conda deactivate
    """
}

workflow alignment_variant_calling {
    take:
    samples
    ref_gz

    main:
    ref = file('ref/GCF_000149955.1_ASM14995v2_genomic.fna')
    fastp(samples)
    Bowtie2(fastp.out.trimmed_reads, ref_gz)
    samTools(Bowtie2.out.sam)
    BCFtools(samTools.out.bam_bai, ref)
    SnpEff(BCFtools.out.vcf)
}