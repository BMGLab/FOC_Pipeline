#!/usr/bin/env nextflow
include { FastQC } from './modules/fastqc.nf'
include { preprocess_assembly } from './workflows/preprocess_assembly.nf'
include { alignment_variant_calling } from './workflows/alignment_variant_calling.nf'
include { QUAST } from './modules/quast.nf'

nextflow.enable.dsl=2

// Initialize the conda environment
// source /home/sercanozturk/miniconda3/etc/profile.d/conda.sh

process inspectReads {
    tag "$id"

    input:
    tuple val(id), path(read1), path(read2)
    path reference

    script:
    """
    echo "Sample ID: $id"
    echo "Read 1: $read1"
    echo "Read 2: $read2"
    """
}

process BLASTn {
    tag "$id"
    publishDir "${params.output_dir}/blastn/${id}", mode: 'copy'

    input:
    val(id)
    path(query)
    // path(blastdb)

    output:
    val(id), emit: sample_id
    path("${query.baseName}_blast_results.txt"), emit: blast_results
    path("${query.baseName}_blast_results.html"), emit: blast_html

    // Coverage of the contigs
    // --qcov_hsp_perc 90.0
    // --perc_identity 98.0
    // -taxids 5507
    script:
    """
    source $params.conda_shell
    conda activate blast
    export BLASTDB=${params.blast_db_dir}
    blastn \
        -query ${query} \
        -db core_nt \
        -taxids 5507 \
        -evalue ${params.blast_evalue} \
        -max_target_seqs ${params.blast_max_target_seqs} \
        -qcov_hsp_perc 90.0 \
        -perc_identity 98.0 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out ${query.baseName}_blast_results.txt

    echo '<html><head><title>BLAST Results</title></head><body>' > ${query.baseName}_blast_results.html
    echo '<h1>BLAST Results</h1>' >> ${query.baseName}_blast_results.html
    echo '<pre>' >> ${query.baseName}_blast_results.html
    cat ${query.baseName}_blast_results.txt | \
        awk 'BEGIN{print "Query\\tSubject\\tIdentity\\tLength\\tMismatches\\tGaps\\tQ.Start\\tQ.End\\tS.Start\\tS.End\\tE-value\\tBit Score"} \
            {print}' | \
        column -t >> ${query.baseName}_blast_results.html
    echo '</pre></body></html>' >> ${query.baseName}_blast_results.html
    conda deactivate
    """
}

process appendContigs {
    tag "$id"
    publishDir "${params.output_dir}/appended_contigs/${id}", mode: 'copy'

    input:
    val id
    path blast_results
    path chr0_fasta
    path ragtag_scaffolds

    output:
    val id, emit: sample_id
    path "${id}_final_scaffolds.fasta", emit: final_scaffolds

    script:
    """
    source $params.conda_shell
    conda activate seqtk
    # Filter the best hit per Chr0 contig by highest bit score
    sort -k12,12nr $blast_results | awk '!seen[\$1]++' > best_hits.tsv

    # Extract only those Chr0 contigs with valid BLAST placements
    awk '{print \$1}' best_hits.tsv > filtered_chr0_ids.txt
    seqtk subseq $chr0_fasta filtered_chr0_ids.txt > filtered_chr0.fasta

    # Remove artificial Chr0 from ragtag scaffolds
    seqkit grep -r -v -p 'Chr0_RagTag' $ragtag_scaffolds > cleaned_scaffolds.fasta

    # Merge the filtered Chr0 contigs with cleaned scaffolds
    cat cleaned_scaffolds.fasta filtered_chr0.fasta > ${id}_final_scaffolds.fasta
    conda deactivate
    """
}

process nucmerMummer {
    tag "$id"
    publishDir "${params.output_dir}/nucmer/${id}", mode: 'copy'

    input:
    val id
    path scaffold
    path reference

    output:
    path "dotplot.png", emit: plot
    path "alignment.delta", emit: delta

    script:
    """
    source $params.conda_shell
    conda activate mummer
    nucmer -p alignment ${scaffold} ${reference}
    mummerplot --png -p dotplot alignment.delta
    conda deactivate
    """
}

process repeatModeler {
    tag "$id"
    publishDir "${params.output_dir}/repeatmodeler/${id}", mode: 'copy'

    input:
    val id
    path scaffold

    output:
    path "query_db-families.fa", emit: repeat_lib
    val id, emit: sample_id

    script:
    """
    source $params.conda_shell
    conda activate repeatmodeler
    BuildDatabase -name query_db ${scaffold}
    RepeatModeler -database query_db -threads 72
    conda deactivate
    """
}

process repeatMasker {
    tag "$id"
    publishDir "${params.output_dir}/repeatmasker/${id}", mode: 'copy'

    input:
    val id
    path scaffold
    path repeat_lib

    output:
    path "${scaffold.baseName}_masked.fasta.gff", emit: masked_gff
    path "${scaffold.baseName}_masked.fasta", emit: masked_fasta
    val id, emit: sample_id

    script:
    """
    source $params.conda_shell
    conda activate repeatmodeler
    RepeatMasker -lib ${repeat_lib} -pa 72 -gff -dir . ${scaffold}
    mv ${scaffold}.masked ${scaffold.baseName}_masked.fasta
    mv ${scaffold}.out.gff ${scaffold.baseName}_masked.fasta.gff
    conda deactivate
    """
}

process summarizeTE {
    tag "$id"
    publishDir "${params.output_dir}/te_summary/${id}", mode: 'copy'

    input:
    val id
    path repeat_lib
    path masked_fasta
    path masked_gff

    output:
    path "${id}_te_family_counts.tsv", emit: family_counts
    path "${id}_te_chrom_density.tsv", emit: chrom_density

    script:
    """
    python - <<'PY'
    import re
    from collections import Counter, defaultdict

    sample_id = "${id}"
    repeat_lib = "${repeat_lib}"
    masked_fasta = "${masked_fasta}"
    masked_gff = "${masked_gff}"

    family_to_class = {}
    with open(repeat_lib, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            left = header.split("(")[0].strip()
            if "#" in left:
                family, klass = left.split("#", 1)
                family_to_class[family.strip()] = klass.strip()

    chrom_lengths = defaultdict(int)
    current = None
    with open(masked_fasta, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].strip().split()[0]
                continue
            if current:
                chrom_lengths[current] += len(line.strip())

    te_types = {"repeat_region", "transposable_element", "inverted_repeat_region"}
    family_counts = Counter()
    class_counts = Counter()
    chrom_te_bp = Counter()

    with open(masked_gff, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, _, _, _, attrs = parts
            if feature not in te_types and source.lower() != "repeatmasker":
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            length = max(0, end_i - start_i + 1)
            chrom_te_bp[chrom] += length

            motif = None
            m = re.search(r"Motif:([^;\\s]+)", attrs)
            if m:
                motif = m.group(1)
            family = motif if motif else "Unknown"
            family_counts[family] += 1

            klass = family_to_class.get(family, "Unknown")
            class_counts[klass] += 1

    with open(f"{sample_id}_te_family_counts.tsv", "w", encoding="utf-8") as out:
        out.write("Sample\\tFamily\\tClass\\tCount\\n")
        for family, count in sorted(family_counts.items(), key=lambda x: (-x[1], x[0])):
            out.write(f"{sample_id}\\t{family}\\t{family_to_class.get(family, 'Unknown')}\\t{count}\\n")

    with open(f"{sample_id}_te_chrom_density.tsv", "w", encoding="utf-8") as out:
        out.write("Sample\\tChromosome\\tChromLengthBp\\tTEBp\\tTEPercent\\n")
        for chrom, clen in sorted(chrom_lengths.items()):
            tebp = chrom_te_bp.get(chrom, 0)
            pct = (100.0 * tebp / clen) if clen else 0.0
            out.write(f"{sample_id}\\t{chrom}\\t{clen}\\t{tebp}\\t{pct:.4f}\\n")
    PY
    """
}

process Liftoff {
    tag "$id"
    publishDir "${params.output_dir}/liftoff/${id}", mode: 'copy'

    input:
    val id
    path assembly
    path reference
    path annotation

    output:
    val id, emit: sample_id
    path "${id}_predicted_genes.gff", emit: genes_gff
    path assembly, emit: scaffold

    script:
    """
    source $params.conda_shell
    conda activate liftoff
    liftoff -g $annotation -o ${id}_predicted_genes.gff \
            -p 28 \
            $assembly $reference
    conda deactivate
    """
}


process SIXBlastp {
    tag "$id"
    publishDir "${params.output_dir}/six_blastp/${id}", mode: 'copy'

    input:
    val id
    path proteins_fasta
    path six_queries

    output:
    path "${id}_six_blastp_raw.tsv", emit: raw
    path "${id}_six_blastp_filtered.tsv", emit: filtered

    script:
    """
    source $params.conda_shell
    conda activate blast

    makeblastdb \
        -in ${proteins_fasta} \
        -dbtype prot \
        -out ${id}_proteome_db

    blastp \
        -query ${six_queries} \
        -db ${id}_proteome_db \
        -evalue ${params.six_blast_evalue} \
        -num_threads ${params.six_blast_threads} \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs" \
        -out ${id}_six_blastp_raw.tsv

    awk -v min_qcov="${params.six_min_query_coverage}" -v max_e="${params.six_blast_evalue}" \
        '(\$11+0) <= (max_e+0) && (\$13+0) >= (min_qcov+0)' \
        ${id}_six_blastp_raw.tsv > ${id}_six_blastp_filtered.tsv

    conda deactivate
    """
}

process dbCAN {
    tag "$id"
    publishDir "${params.output_dir}/dbcan/${id}", mode: 'copy'

    input:
    val id
    path proteins

    output:
    val id, emit: sample_id
    path "out", emit: cazymes
    path "${id}_dbcan_hmmer_filtered.tsv", emit: hmmer_filtered
    path "${id}_dbcan_protein_summary.tsv", emit: protein_summary

    script:
    """
    source $params.conda_shell
    conda activate dbcan

    run_dbcan CAZyme_annotation \
        --mode protein \
        --input_raw_data ${proteins} \
        --output_dir out \
        --db_dir $params.dbcan_db_dir \
        --methods hmm \
        --hmm_eval ${params.dbcan_hmm_dome} \
        --threads ${params.dbcan_threads}

    HMM_OUT=""
    for candidate in out/hmmer.out out/hmm.out out/hmmscan.out; do
        if [[ -s "\$candidate" ]]; then
            HMM_OUT="\$candidate"
            break
        fi
    done
    [[ -n "\$HMM_OUT" ]] || { echo "[dbCAN] HMM output not found in out/"; exit 1; }

    awk -v OFS='\\t' -v thr='${params.dbcan_hmm_dome}' '
        BEGIN { print "ProteinID","CAZyFamily","domE" }
        !/^#/ && NF>=5 && \$5 ~ /^([0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)\$/ {
            family=\$1; protein=\$3; dome=\$5+0;
            if (dome <= thr+0) print protein, family, dome;
        }
    ' "\$HMM_OUT" > ${id}_dbcan_hmmer_filtered.tsv

    awk -v OFS='\\t' '
        BEGIN { print "ProteinID","CAZyFamilyCount","CAZyFamilies" }
        NR>1 {
            protein=\$1; family=\$2;
            key=protein FS family;
            if (!seen[key]++) {
                count[protein]++;
                fams[protein] = (fams[protein] ? fams[protein]","family : family);
            }
        }
        END {
            for (p in count) print p, count[p], fams[p];
        }
    ' ${id}_dbcan_hmmer_filtered.tsv | sort -k1,1 > ${id}_dbcan_protein_summary.tsv

    conda deactivate
    """
}

process extractProteins {
    tag "$id"

    input:
    val id
    path assembly
    path gff

    output:
    val id, emit: sample_id
    path "${id}_predicted_proteins.faa", emit: proteins

    script:
    """
    source $params.conda_shell
    conda activate agat

    agat_sp_fix_cds_phases.pl \
        --gff $gff \
        --fasta ${assembly} \
        -o ${id}_predicted_genes_fixed.gff

    agat_sp_extract_sequences.pl \
        --gff ${id}_predicted_genes_fixed.gff \
        --fasta ${assembly} \
        --type CDS \
        --protein \
        -o ${id}_predicted_proteins.faa

    conda deactivate
    """
}

process TargetP {
    tag "$id"
    publishDir "${params.output_dir}/targetp/${id}", mode: 'copy'

    input:
    val id
    path protein_fasta

    output:
    path("${id}_targetp.gff3"), emit: targetp_gff

    script:
    """
    source $params.conda_shell
    $params.targetp \
        -fasta $protein_fasta \
        -format short \
        -org non-pl \
        -gff3 \
        -mature \
        -prefix ${id}_targetp
    """
}

process Signalp {
    tag "$id"
    publishDir "${params.output_dir}/signalp/${id}", mode: 'copy'

    input:
    val id
    path protein_fasta

    output:
    path("output/*"), emit: signalp_output
    path("output/prediction_results.txt"), emit: prediction_results

    script:
    """
    source $params.conda_shell
    conda activate signalp6
    signalp6 \
        --fastafile $protein_fasta \
        --output_dir output \
        --format all \
        --organism euk \
        --mode slow \
        --model_dir ${params.home}/signalp6_slow_sequential/signalp-6-package/models/
    conda deactivate
    """
}

process WoLFPSort {
    tag "$id"
    publishDir "${params.output_dir}/wolfpsort/${id}", mode: 'copy'

    input:
    val id
    path protein_fasta

    output:
    path("${id}_wolfpsort.txt"), emit: wolfpsort_out

    script:
    """
    ${params.home}/WoLFPSort/bin/runWolfPsortSummary fungi < $protein_fasta > ${id}_wolfpsort.txt
    """
}

workflow {
    if (!params.enable_chr0_blast_append) { println "INFO: Using RagTag scaffold directly (manuscript-aligned; Chr0 BLAST append disabled)\n" }
    if (!params.run_mcclintock) { println "INFO: Skipping McClintock (not in manuscript core TE workflow)\n" }
    println "INFO: Option A enabled: Liftoff annotations are the source of truth for downstream protein/functional analyses\n"
    channel.fromPath(params.samplesheet_path)
        .splitCsv(header: true)
        .map {row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
        .set { samples }
    six_queries = file(params.six_queries_fasta)
    if (!six_queries.exists()) {
        throw new IllegalArgumentException("SIX query FASTA not found at --six_queries_fasta: ${params.six_queries_fasta}")
    }

    reference = file(params.reference_genome)
    reference_fna = file("ref/GCF_000149955.1_ASM14995v2_genomic.fna")
    annotation = file("ref/GCF_000149955.1_ASM14995v2_genomic.gff")

    inspectReads(samples, reference)
    FastQC(samples.merge(channel.value(".")))

    preprocess_assembly(samples, reference)
    alignment_variant_calling(samples, reference)

    def assembly_sample_ids
    def assembly_scaffolds
    if (params.enable_chr0_blast_append) {
        BLASTn(preprocess_assembly.out.sample_id, preprocess_assembly.out.chr0_contigs)
        appendContigs(BLASTn.out.sample_id, BLASTn.out.blast_results, preprocess_assembly.out.chr0_contigs, preprocess_assembly.out.scaffold)
        assembly_sample_ids = appendContigs.out.sample_id
        assembly_scaffolds = appendContigs.out.final_scaffolds
    } else {
        assembly_sample_ids = preprocess_assembly.out.sample_id
        assembly_scaffolds = preprocess_assembly.out.scaffold
    }

    QUAST(assembly_sample_ids, assembly_scaffolds, channel.value("."))
    Liftoff(assembly_sample_ids, assembly_scaffolds, reference_fna, annotation)
    extractProteins(Liftoff.out.sample_id, Liftoff.out.scaffold, Liftoff.out.genes_gff)
    SIXBlastp(extractProteins.out.sample_id, extractProteins.out.proteins, six_queries)
    dbCAN(extractProteins.out.sample_id, extractProteins.out.proteins)

    nucmerMummer(assembly_sample_ids, assembly_scaffolds, reference_fna)
    repeatModeler(assembly_sample_ids, assembly_scaffolds)
    repeatMasker(repeatModeler.out.sample_id, assembly_scaffolds, repeatModeler.out.repeat_lib)
    summarizeTE(repeatMasker.out.sample_id, repeatModeler.out.repeat_lib, repeatMasker.out.masked_fasta, repeatMasker.out.masked_gff)

    TargetP(extractProteins.out.sample_id, extractProteins.out.proteins)
    Signalp(extractProteins.out.sample_id, extractProteins.out.proteins)
    WoLFPSort(extractProteins.out.sample_id, extractProteins.out.proteins)
}
