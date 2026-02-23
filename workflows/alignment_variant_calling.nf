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

process GroupSpecificVariants {
    publishDir "${params.output_dir}/group_specific_variants", mode: 'copy'

    input:
    path callable_vcfs
    path group_map

    output:
    path "merged_callable_variants.vcf", emit: merged_vcf
    path "group_specific_summary.tsv", emit: summary
    path "group_specific/*.vcf", emit: group_vcfs

    script:
    """
    source $params.conda_shell
    conda activate bcftools

    renamed_vcfs=()
    for vcf in ${callable_vcfs}; do
        id=\$(basename "\$vcf" | sed 's/_callable_filtered_variants\\.vcf\$//')
        echo "\$id" > "\${id}.sample_name.txt"
        bcftools reheader -s "\${id}.sample_name.txt" -o "\${id}.renamed.vcf" "\$vcf"
        renamed_vcfs+=("\${id}.renamed.vcf")
    done

    bcftools merge -Ov -o merged_callable_variants.vcf "\${renamed_vcfs[@]}"
    conda deactivate

    python - <<'PY'
    from collections import defaultdict
    from pathlib import Path

    group_map_path = Path("${group_map}")
    merged_vcf_path = Path("merged_callable_variants.vcf")
    out_dir = Path("group_specific")
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_to_group = {}
    with group_map_path.open("r", encoding="utf-8") as fh:
        for i, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\\t")
            if len(cols) < 2:
                raise SystemExit(f"Invalid group map at line {i}: expected 'sample_id<TAB>group'")
            sample_to_group[cols[0]] = cols[1]

    meta = []
    header = None
    samples = []
    with merged_vcf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("##"):
                meta.append(line)
                continue
            if line.startswith("#CHROM"):
                header = line
                samples = line.rstrip("\\n").split("\\t")[9:]
                break

    if not samples:
        raise SystemExit("No samples found in merged VCF.")

    missing = [s for s in samples if s not in sample_to_group]
    if missing:
        raise SystemExit(f"Samples missing in group map: {', '.join(missing)}")

    groups = sorted(set(sample_to_group[s] for s in samples))
    if len(groups) < 2:
        raise SystemExit("Need at least two groups for strict group-specific variant calling.")

    group_samples = {g: [s for s in samples if sample_to_group[s] == g] for g in groups}
    counts = defaultdict(int)

    handles = {}
    for g in groups:
        out_path = out_dir / f"{g}_strict_specific.vcf"
        h = out_path.open("w", encoding="utf-8")
        for m in meta:
            h.write(m)
        h.write(header)
        handles[g] = h

    def parse_gt(sample_field, gt_idx):
        parts = sample_field.split(":")
        if gt_idx >= len(parts):
            return "./."
        return parts[gt_idx]

    def is_callable(gt):
        return gt not in {".", "./.", ".|."}

    def has_alt(gt):
        if not is_callable(gt):
            return False
        alleles = gt.replace("|", "/").split("/")
        return all(a != "." for a in alleles) and any(a != "0" for a in alleles)

    def is_ref(gt):
        if not is_callable(gt):
            return False
        alleles = gt.replace("|", "/").split("/")
        return all(a == "0" for a in alleles)

    with merged_vcf_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 10:
                continue

            fmt_keys = cols[8].split(":")
            if "GT" not in fmt_keys:
                continue
            gt_idx = fmt_keys.index("GT")
            sample_fields = cols[9:]
            gt_by_sample = {
                sample: parse_gt(field, gt_idx)
                for sample, field in zip(samples, sample_fields)
            }

            matched_group = None
            for g in groups:
                in_group = group_samples[g]
                out_group = [s for s in samples if s not in in_group]
                if not in_group or not out_group:
                    continue
                if all(has_alt(gt_by_sample[s]) for s in in_group) and \
                   all(is_ref(gt_by_sample[s]) for s in out_group):
                    matched_group = g
                    break

            if matched_group:
                handles[matched_group].write(line)
                counts[matched_group] += 1

    for h in handles.values():
        h.close()

    with Path("group_specific_summary.tsv").open("w", encoding="utf-8") as out:
        out.write("Group\\tNumSamples\\tStrictSpecificVariantCount\\n")
        for g in groups:
            out.write(f"{g}\\t{len(group_samples[g])}\\t{counts[g]}\\n")
    PY
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
    group_map = file(params.avc_group_map)
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
    if (params.avc_run_group_specific) {
        if (!group_map.exists()) {
            throw new IllegalArgumentException("Group map not found at --avc_group_map: ${params.avc_group_map}")
        }
        GroupSpecificVariants(RestrictVariantsToCallable.out.vcf.map { id, vcf -> vcf }.collect(), group_map)
    }
    SnpEff(RestrictVariantsToCallable.out.vcf)
}
