include { fastp } from '../modules/fastp.nf'
include { FastQC } from '../modules/fastqc.nf'
include { QUAST } from '../modules/quast.nf'
include { LAST } from '../modules/last.nf'

process megahitAssembly {
    tag "$id"
    publishDir "${params.output_dir}/assembly/${id}", mode: 'copy'

    input:
    tuple val(id), val(read1), val(read2)

    output:
    val id, emit: sample_id
    path "megahit/final.contigs.fa", emit: contigs

    script:
    """
    source ${params.conda_shell}
    conda activate megahit
    megahit -1 '${read1}' -2 '${read2}' \
            --out-dir 'megahit' \
            --presets meta-sensitive \
            --num-cpu-threads 84
    conda deactivate
    """
}

process ragtagCorrect {
    tag "$id"
    publishDir "${params.output_dir}/ragtag/${id}/correction", mode: 'copy'

    input:
    val id
    path contigs
    path reference

    output:
    val id, emit: sample_id
    path "ragtag_correction/ragtag.correct.fasta", emit: corrected_contigs
    path reference

    script:
    """
    source $params.conda_shell
    conda activate ragtag
    ragtag.py correct ${reference} ${contigs} -o ragtag_correction
    conda deactivate
    """
}

process ragtagScaffold {
    tag "$id"
    publishDir "${params.output_dir}/ragtag/${id}/scaffold", mode: 'copy'

    input:
    val id
    path corrected_contigs
    path reference

    output:
    val id, emit: sample_id
    path "minimapRagTag/ragtag.scaffold.fasta", emit: scaffold_fasta
    path "minimapRagTag/ragtag.scaffold_unplaced.fasta", emit: unplaced
    path "minimapRagTag/ragtag.scaffold.agp", emit: scaffold_agp
    path corrected_contigs, emit: corrected_contigs

    script:
    """
    source $params.conda_shell
    conda activate ragtag
    ragtag.py scaffold_custom ${reference} ${corrected_contigs} \
        -f 50 --remove-small -C \
        -o minimapRagTag
    echo ${params.project_root}
    python ${params.project_root}/utilities/agp_to_fasta.py \
    minimapRagTag/ragtag.scaffold_unplaced.agp \
    ${corrected_contigs} \
    minimapRagTag/ragtag.scaffold_unplaced.fasta
    conda deactivate
    """
}

process extractChr0Contigs {
    tag "$id"
    publishDir "${params.output_dir}/chr0_contigs/${id}", mode: 'copy'

    input:
    val id
    path agp_file
    path corrected_fasta

    output:
    path 'chr0_contigs.fasta', emit: chr0_contigs
    path agp_file
    path corrected_fasta
    path 'chr0_contig_ids.txt'

    script:
    """
    source $params.conda_shell
    conda activate ragtag
    grep "Chr0_RagTag" "${agp_file}" | awk '\$5 != "U" {print \$6}' | sort | uniq > "chr0_contig_ids.txt"
    CHR0_COUNT=\$(wc -l < "chr0_contig_ids.txt")
    echo "Found \$CHR0_COUNT Chr0 contigs in AGP file"

    if [[ \$CHR0_COUNT -eq 0 ]]; then
        echo "No Chr0 contigs found. Please check your RagTag output."
        exit 1
    fi

    seqtk subseq "${corrected_fasta}" "chr0_contig_ids.txt" > "chr0_contigs.fasta"
    conda deactivate
    """
}

process bwaCoverageEstimation {
    tag "$id"
    publishDir "${params.output_dir}/coverage/${id}", mode: 'copy'

    input:
    tuple val(id), path(read1), path(read2)
    path reference

    output:
    val id, emit: sample_id
    path "${id}.dedup.bam", emit: bam
    path "${id}.dedup.bam.bai", emit: bai
    path "${id}.depth.tsv", emit: depth
    path "${id}.coverage_metrics.tsv", emit: metrics

    script:
    """
    source $params.conda_shell
    conda activate bwa

    cp ${reference} reference.fa
    bwa index reference.fa
    bwa mem -t ${params.coverage_threads} reference.fa ${read1} ${read2} > ${id}.sam

    conda deactivate
    conda activate samtools

    samtools view -@ ${params.coverage_threads} -b ${id}.sam | \
        samtools sort -@ ${params.coverage_threads} -n -o ${id}.namesort.bam
    samtools fixmate -@ ${params.coverage_threads} -m ${id}.namesort.bam ${id}.fixmate.bam
    samtools sort -@ ${params.coverage_threads} -o ${id}.sorted.bam ${id}.fixmate.bam
    samtools markdup -@ ${params.coverage_threads} -s ${id}.sorted.bam ${id}.dedup.bam
    samtools index -@ ${params.coverage_threads} ${id}.dedup.bam

    samtools depth -a -@ ${params.coverage_threads} ${id}.dedup.bam > ${id}.depth.tsv

    awk -v sample="${id}" '{
        tot++;
        sum += \$3;
        if (\$3 > 0) cov++;
    } END {
        avg = (tot > 0 ? sum / tot : 0);
        breadth = (tot > 0 ? 100 * cov / tot : 0);
        printf "%s\\t%.2f\\t%.2f\\n", sample, avg, breadth;
    }' ${id}.depth.tsv > ${id}.coverage_metrics.tsv

    rm -f ${id}.sam ${id}.namesort.bam ${id}.fixmate.bam ${id}.sorted.bam \
          reference.fa reference.fa.amb reference.fa.ann reference.fa.bwt reference.fa.pac reference.fa.sa
    conda deactivate
    """
}

process summarizeCoverageMetrics {
    publishDir "${params.output_dir}/coverage", mode: 'copy'

    input:
    path metrics_files

    output:
    path "coverage_summary.tsv", emit: summary

    script:
    """
    echo -e "Sample\\tAvgDepthX\\tBreadthPct" > coverage_summary.tsv
    for f in ${metrics_files}; do
        cat "\$f" >> coverage_summary.tsv
    done
    """
}

workflow preprocess_assembly {
    take:
    samples
    reference

    main:
    def output_dir = Channel.value("${params.pa_wf_output}")
    samples.merge(output_dir).set{ samples_ch }
    fastp(samples_ch)
    FastQC(samples_ch)
    megahit_out = megahitAssembly(fastp.out.trimmed_reads)
    QUAST(megahit_out.sample_id, megahit_out.contigs, output_dir)
    ragtag_out = ragtagCorrect(megahit_out.sample_id, megahit_out.contigs, reference) | ragtagScaffold
    LAST(ragtag_out.sample_id, ragtag_out.scaffold_fasta, reference)
    bwaCoverageEstimation(fastp.out.trimmed_reads, reference)
    summarizeCoverageMetrics(bwaCoverageEstimation.out.metrics.collect())
    def chr0_contigs_out = Channel.empty()
    if (params.enable_chr0_blast_append) {
        extractChr0Contigs(ragtag_out.sample_id, ragtag_out.scaffold_agp, ragtag_out.corrected_contigs)
        chr0_contigs_out = extractChr0Contigs.out.chr0_contigs
    }

    emit:
    sample_id = ragtag_out.sample_id
    chr0_contigs = chr0_contigs_out
    scaffold = ragtag_out.scaffold_fasta
    coverage_metrics = bwaCoverageEstimation.out.metrics
    coverage_summary = summarizeCoverageMetrics.out.summary
}
