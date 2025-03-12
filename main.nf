#!/usr/bin/env nextflow
include { fastqc } from './modules/fastqc.nf'
include { preprocess_assembly } from './workflows/preprocess_assembly.nf'
include { alignment_variant_calling } from './workflows/alignment_variant_calling.nf'

nextflow.enable.dsl=2

// Initialize the conda environment
// source /home/sercanozturk/miniconda3/etc/profile.d/conda.sh

process inspectReads {
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

process makeBLASTDB {
    publishDir "${params.output_dir}/blast_db", mode: 'copy'
    
    input:
    path reference
    
    output:
    path "blastdb*", emit: blast_db
    
    script:
    """
    source $params.conda_shell
    conda activate blast
    makeblastdb \
        -in ${reference} \
        -dbtype nucl \
        -out blastdb \
        -parse_seqids
    conda deactivate
    """
}

process BLASTSearch {
    tag "$id"
    publishDir "${params.output_dir}/blast/${id}", mode: 'copy'

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

    input:
    val id
    path scaffold
    path unplaced
    path blast

    output:
    val(id), emit: sample_id
    path("final_scaffolds.fasta"), emit: final_scaffolds

    script:
    """
    grep -E "^>" ${blast} | cut -f1 -d' ' > fusarium_contigs.txt
    seqkit grep -f fusarium_contigs.txt ${unplaced} > fusarium_contigs.fasta
    cat ${scaffold} fusarium_contigs.fasta > final_scaffolds.fasta
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
    path "RepeatModelerLib.fa", emit: repeat_lib
    val id, emit: sample_id

    script:
    """
    source $params.conda_shell
    conda activate repeatmodeler
    BuildDatabase -name query_db ${scaffold}
    RepeatModeler -database query_db -threads ${task.cpus}
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

    script:
    """
    source $params.conda_shell
    conda activate repeatmodeler
    RepeatMasker -lib ${repeat_lib} -gff -dir . ${scaffold}
    mv ${scaffold}.masked ${scaffold.baseName}_masked.fasta
    mv ${scaffold}.gff ${scaffold.baseName}_masked.fasta.gff
    conda deactivate
    """
}

process mcclintock {
    tag "$id"
    publishDir "${params.output_dir}/mcclintock/${id}", mode: 'copy'

    input:
    val id
    path masked_gff
    path scaffold

    output:
    path "mcclintock_results/*", emit: results

    script:
    """
    source $params.conda_shell
    conda activate mcclintock
    mcclintock --gff ${masked_gff} \
               --fasta ${scaffold} \
               --output_dir mcclintock_results
    conda deactivate
    """
}

workflow {
    Channel.fromPath(params.samplesheet_path).splitCsv(header: true).map {row -> 
            tuple(row.sample_id, file(row.read1), row.read2)}.set { samples }

    reference = file(params.reference_genome)
    reference_fna = file("ref/GCF_000149955.1_ASM14995v2_genomic.fna")
    inspectReads(samples, reference)
    fastqc(samples)

    preprocess_assembly(samples, reference)
    alignment_variant_calling(samples, reference)

    // makeBLASTDB(reference_fna)   -- No need to make when using core_nt
    BLASTSearch(preprocess_assembly.out.sample_id, preprocess_assembly.out.chr0_contigs)

    nucmerMummer(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold, reference_fna)
    repeatModeler(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold)
    repeatMasker(repeatModeler.out.sample_id, preprocess_assembly.out.scaffold, repeatModeler.out.repeat_lib)
    mcclintock(BLASTSearch.out.sample_id, repeatMasker.out.masked_gff, BLASTSearch.out.blast_results)
}