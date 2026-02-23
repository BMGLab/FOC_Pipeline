include { fastp } from '../modules/fastp.nf'
include { samTools } from '../modules/samtools.nf'

process Bowtie2 {
    tag "$id"
    publishDir "${params.output_dir}/alignment/bowtie2/${id}", mode: 'copy'

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
            2> ${id}_bowtie2.log \
            --threads 28
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
    bcftools mpileup \
        -f ${reference} \
        --max-depth ${params.avc_mpileup_max_depth} \
        -q ${params.avc_min_mapq} \
        -Q ${params.avc_min_baseq} \
        ${bam} -Ou | \
    bcftools call \
        --ploidy ${params.avc_ploidy} \
        -mv \
        -Ov \
        -o ${id}_variants.bcf

    bcftools filter \
        -s LowQual \
        -e 'QUAL<${params.avc_variant_min_qual} || INFO/DP<${params.avc_variant_min_dp}' \
        -o ${id}_filtered_variants.bcf \
        ${id}_variants.bcf

    bcftools view ${id}_filtered_variants.bcf -Ov -o ${id}_filtered_variants.vcf
    conda deactivate
    """
}

process CallableRegions {
    tag "$id"
    publishDir "${params.output_dir}/callable/${id}", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)
    path reference

    output:
    tuple val(id), path("${id}_callable.bed"), emit: callable_bed
    path "${id}.mosdepth.summary.txt", emit: mosdepth_summary

    script:
    """
    source $params.conda_shell
    conda activate mosdepth
    mosdepth \
        --threads ${params.avc_threads} \
        --mapq ${params.avc_min_mapq} \
        ${id} \
        ${bam}
    conda deactivate

    source $params.conda_shell
    conda activate samtools
    samtools depth \
        -aa \
        -q ${params.avc_min_mapq} \
        -Q ${params.avc_min_baseq} \
        ${bam} | \
    awk -v min_cov=${params.avc_min_callable_depth} 'BEGIN{OFS="\\t"}
        \$3 >= min_cov {
            chr=\$1; pos=\$2;
            if (current_chr != chr || pos != prev_pos + 1) {
                if (current_chr != "") print current_chr, start_pos - 1, prev_pos;
                current_chr = chr;
                start_pos = pos;
            }
            prev_pos = pos;
        }
        END {
            if (current_chr != "") print current_chr, start_pos - 1, prev_pos;
        }' > ${id}_callable.bed
    conda deactivate
    """
}

process ConsensusCallableRegions {
    publishDir "${params.output_dir}/callable", mode: 'copy'

    input:
    path callable_beds
    path reference

    output:
    path "consensus_callable.bed", emit: consensus_bed
    path "consensus_callable_stats.tsv", emit: consensus_stats

    script:
    """
    source $params.conda_shell
    conda activate samtools
    samtools faidx ${reference}
    conda deactivate

    python - <<'PY'
    import math
    from collections import defaultdict

    ref_fai = "${reference}.fai"
    bed_files = "${callable_beds}".split()
    frac = float("${params.avc_consensus_callable_fraction}")

    chrom_lengths = {}
    with open(ref_fai, "r", encoding="utf-8") as fh:
        for line in fh:
            cols = line.rstrip("\\n").split("\\t")
            chrom_lengths[cols[0]] = int(cols[1])

    n = len(bed_files)
    if n == 0:
        raise SystemExit("No callable BED files were provided.")
    k = max(1, math.ceil(n * frac))

    events = defaultdict(list)
    for bed in bed_files:
        with open(bed, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                chrom, s, e = line.rstrip("\\n").split("\\t")[:3]
                start = int(s)
                end = int(e)
                if chrom not in chrom_lengths:
                    continue
                events[chrom].append((start, 1))
                events[chrom].append((end, -1))

    total_callable = 0
    total_genome = sum(chrom_lengths.values())
    with open("consensus_callable.bed", "w", encoding="utf-8") as out_bed:
        for chrom, clen in chrom_lengths.items():
            chrom_events = sorted(events.get(chrom, []))
            if not chrom_events:
                continue
            cov = 0
            prev = 0
            for pos, delta in chrom_events:
                if pos > prev and cov >= k:
                    out_bed.write(f"{chrom}\\t{prev}\\t{pos}\\n")
                    total_callable += (pos - prev)
                cov += delta
                prev = pos
            if prev < clen and cov >= k:
                out_bed.write(f"{chrom}\\t{prev}\\t{clen}\\n")
                total_callable += (clen - prev)

    pct = (100.0 * total_callable / total_genome) if total_genome else 0.0
    with open("consensus_callable_stats.tsv", "w", encoding="utf-8") as out_stats:
        out_stats.write("NumSamples\\tKThreshold\\tFraction\\tCallableBp\\tGenomeBp\\tCallablePct\\n")
        out_stats.write(f"{n}\\t{k}\\t{frac:.4f}\\t{total_callable}\\t{total_genome}\\t{pct:.4f}\\n")
    PY
    """
}

process RestrictVariantsToCallable {
    tag "$id"
    publishDir "${params.output_dir}/bcftools_callable/${id}", mode: 'copy'

    input:
    tuple val(id), path(vcf)
    path consensus_bed

    output:
    tuple val(id), path("${id}_callable_filtered_variants.vcf"), emit: vcf

    script:
    """
    source $params.conda_shell
    conda activate bcftools
    bcftools view -R ${consensus_bed} ${vcf} -Ov -o ${id}_callable_filtered_variants.vcf
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

    conda activate bcftools
    bcftools annotate \
        --rename-chrs ${params.chr_rename_map} \
        --output ${id}_renamed.vcf \
        ${vcf}
    conda deactivate

    conda activate snpeff
    export SNPEFF_JAR=${params.snpeff_jar}
    java -jar \$SNPEFF_JAR databases | grep ${params.snpeff_db}
    java -Xmx8g -jar \$SNPEFF_JAR \
        -v ${params.snpeff_db} \
        -stats ${id}_snpEff_summary.html \
        -csvStats snpEff_summary.csv \
        ${id}_renamed.vcf > ${id}_annotated.vcf
    conda deactivate

    grep -v "^#" ${id}_annotated.vcf | cut -f 8 | grep -o 'ANN=.*' | \
    sed 's/ANN=//g' | tr ',' '\\n' | cut -f 4 -d '|' | sort | uniq -c > snpEff_genes.txt
    """
}

workflow alignment_variant_calling {
    take:
    samples
    ref_gz

    main:
    ref = file('ref/GCF_000149955.1_ASM14995v2_genomic.fna')
    def output_dir = channel.value("${params.avc_wf_output}")
    samples.merge(output_dir).set { samples_ch }
    fastp(samples_ch)
    Bowtie2(fastp.out.trimmed_reads, ref_gz)
    Bowtie2.out.sam.merge(output_dir).set { bowtie2_ch }
    samTools(bowtie2_ch)
    CallableRegions(samTools.out.bam_bai, ref)
    ConsensusCallableRegions(CallableRegions.out.callable_bed.collect(), ref)
    BCFtools(samTools.out.bam_bai, ref)
    RestrictVariantsToCallable(BCFtools.out.vcf, ConsensusCallableRegions.out.consensus_bed)
    SnpEff(RestrictVariantsToCallable.out.vcf)
}
