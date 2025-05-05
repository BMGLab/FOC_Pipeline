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

// Need to bypass "Is out folder empty?" check in mcclintock.py
process McClintock {
    tag "$id"
    publishDir "${params.output_dir}/mcclintock/${id}", mode: 'copy'

    input:
    val id
    path reference
    path consensus
    path read1
    path read2
    path locations

    output:
    path "out/*", emit: results

    script:
    def read1_name_fq = read1.getBaseName()
    def read1_name = read1_name_fq.substring(0, read1_name_fq.lastIndexOf('.'))
    """
    source $params.conda_shell
    conda activate mcclintock

    python - <<EOF
    import re
    def clean_fasta_headers(input_file):
        output_file = input_file + "_cleaned"
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    clean_header = re.sub(r"[^a-zA-Z0-9_.#>]", "_", line.strip())
                    outfile.write(clean_header + "\\n")
                else:
                    outfile.write(line)
        return output_file
    cleaned_consensus = clean_fasta_headers("${consensus}")
    EOF

    mkdir -p out/${read1_name}/results/relocate2/unfiltered/repeat/results
    python $params.mcclintock \
        --reference ${reference} \
        --consensus ${consensus}_cleaned \
        --first ${read1} \
        --second ${read2} \
        --proc 28 \
        --out out \
        --methods ngs_te_mapper2,relocate2,tebreak \
        --keep_intermediate all
    conda deactivate
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

process AUGUSTUS {
    tag "$id"
    publishDir "${params.output_dir}/augustus/${id}", mode: 'copy'

    input:
    val id
    path assembly

    output:
    val id, emit: sample_id
    path "${id}_augustus_predictions.gff", emit: genes_gff
    path assembly, emit: scaffold
    path "${id}_augustus_predictions.faa", emit: proteins_faa

    script:
    """
    source $params.conda_shell
    conda activate augustus

    augustus \
        --species=fol \
        --gff3=on \
        --outfile=${id}_augustus_predictions.gff \
        $assembly

    getAnnoFasta.pl ${id}_augustus_predictions.gff
    mv ${id}_augustus_predictions.aa ${id}_augustus_predictions.faa

    conda deactivate
    """
}

process antiSMASH {
    tag "$id"
    publishDir "${params.output_dir}/antismash/${id}", mode: 'copy'

    input:
    val id
    path assembly
    path gff3

    output:
    path "${id}_antismash_results/", emit: antismash

    script:
    """
    source $params.conda_shell
    conda activate agat
    agat_sp_manage_IDs.pl --gff ${gff3} -o ${id}_fixed_genes.gff
    agat_sp_fix_cds_phases.pl --gff ${id}_fixed_genes.gff --fasta ${assembly} -o ${id}_fixed_genes_phased.gff
    agat_sp_fix_features_locations_duplicated.pl --gff ${id}_fixed_genes_phased.gff -o ${id}_fixed_genes_phased_deduped.gff
    agat_sp_keep_longest_isoform.pl -g ${id}_fixed_genes_phased_deduped.gff -o no_redundant_isoforms.gff
    conda deactivate
    conda activate antismash
    antismash --taxon fungi \
                --cpus 28 \
                --output-dir ${id}_antismash_results \
                --fullhmmer --clusterhmmer --cc-mibig --cb-general --cb-knownclusters \
                --genefinding-gff3 no_redundant_isoforms.gff \
                ${assembly}
    conda deactivate
    """
}

process BLASTp {
    tag "$id"
    publishDir "${params.output_dir}/blastp/${id}", mode: 'copy'

    input:
    val id
    path proteins_fasta

    output:
    path "${id}_blastp_merops.txt", emit: blastp_merops

    script:
    def merops_db = file("${params.merops_db_dir}")
    """
    source $params.conda_shell
    conda activate blast
    blastp -query $proteins_fasta \
           -db ${merops_db} \
           -out ${id}_blastp_merops.txt \
           -evalue 1e-5 \
           -num_threads 28 \
           -max_target_seqs 10 \
           -outfmt 6
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

    script:
    """
    source $params.conda_shell
    conda activate dbcan

    run_dbcan CAZyme_annotation \
        --mode protein \
        --input_raw_data ${proteins} \
        --output_dir out \
        --db_dir $params.dbcan_db_dir \
        --methods diamond,hmm,dbCANsub \
        --threads 4

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

process DeepTMHMM {
    tag "$id"
    publishDir "${params.output_dir}/deeptmhmm/${id}", mode: 'copy'

    input:
    val id
    path protein_fasta

    output:
    path("${id}_tmhmm.out"), emit: tmhmm_out

    script:
    """
    source $params.conda_shell
    conda activate deeptmhmm

    biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta ${protein_fasta} > ${id}_tmhmm.out

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
    path("output/*"), emit: signalp_out

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
        --torch_num_threads 28 \
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
    if (params.skip_deeptmhmm) { println "INFO: Skipping DeepTMHMM\n" }
    Channel.fromPath(params.samplesheet_path)
        .splitCsv(header: true)
        .map {row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
        .set { samples }

    reference = file(params.reference_genome)
    reference_fna = file("ref/GCF_000149955.1_ASM14995v2_genomic.fna")
    annotation = file("ref/GCF_000149955.1_ASM14995v2_genomic.gff")

    inspectReads(samples, reference)
    FastQC(samples.merge(Channel.value(".")))

    preprocess_assembly(samples, reference)
    alignment_variant_calling(samples, reference)

    BLASTn(preprocess_assembly.out.sample_id, preprocess_assembly.out.chr0_contigs)
    appendContigs(BLASTn.out.sample_id, BLASTn.out.blast_results, preprocess_assembly.out.chr0_contigs, preprocess_assembly.out.scaffold)
    QUAST(appendContigs.out.sample_id, appendContigs.out.final_scaffolds, Channel.value("."))
    Liftoff(appendContigs.out.sample_id, appendContigs.out.final_scaffolds, reference_fna, annotation)
    AUGUSTUS(appendContigs.out.sample_id, appendContigs.out.final_scaffolds)
    dbCAN(AUGUSTUS.out.sample_id, AUGUSTUS.out.proteins_faa)
    BLASTp(AUGUSTUS.out.sample_id, AUGUSTUS.out.proteins_faa)
    antiSMASH(AUGUSTUS.out.sample_id, AUGUSTUS.out.scaffold, AUGUSTUS.out.genes_gff)

    nucmerMummer(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold, reference_fna)
    repeatModeler(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold)
    repeatMasker(repeatModeler.out.sample_id, preprocess_assembly.out.scaffold, repeatModeler.out.repeat_lib)
    read1 = samples.map { it[1] }
    read2 = samples.map { it[2] }
    McClintock(repeatMasker.out.sample_id, repeatMasker.out.masked_fasta, repeatModeler.out.repeat_lib, read1, read2, repeatMasker.out.masked_gff)

    extractProteins(Liftoff.out.sample_id, Liftoff.out.scaffold, Liftoff.out.genes_gff)
    
    if(!params.skip_deeptmhmm) { DeepTMHMM(extractProteins.out.sample_id, extractProteins.out.proteins) }
    TargetP(extractProteins.out.sample_id, extractProteins.out.proteins)
    Signalp(extractProteins.out.sample_id, extractProteins.out.proteins)
    WoLFPSort(extractProteins.out.sample_id, extractProteins.out.proteins)
}
