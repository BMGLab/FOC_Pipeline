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

process mcclintock {
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
    path "mcclintock_results/*", emit: results

    script:
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

    python /home/sercanozturk/Fol_Nextflow_Pipeline/modules/McClintock/mcclintock.py \
        --reference ${reference} \
        --consensus ${consensus}_cleaned \
        --first ${read1} \
        --second ${read2} \
        --out mcclintock_results \
        --proc 12
    conda deactivate
    """
}

process liftoff {
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

process antismash {
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
    conda activate antismash
    antismash --taxon fungi \
                --cpus 28 \
                --output-dir ${id}_antismash_results \
                --fullhmmer --clusterhmmer --cc-mibig --cb-general --cb-knownclusters \
                --genefinding-gff3 ${gff3} \
                ${assembly}
    conda deactivate
    """
}

process blastp {
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

process dbscan {
    tag "$id"
    publishDir "${params.output_dir}/dbscan/${id}", mode: 'copy'

    input:
    val id
    path gff

    output:
    val id, emit: sample_id
    path "clustered_genes.tsv", emit: clustered_genes

    script:
    """
    source $params.conda_shell
    conda activate dbscan
    python ${params.project_root}/modules/dbscan.py \
        --input ${gff} \
        --output clustered_genes.tsv \
        --eps 1000 \
        --min_samples 3
    conda deactivate
    """
}

workflow {
    Channel.fromPath(params.samplesheet_path)
        .splitCsv(header: true)
        .map {row -> tuple(row.sample_id, file(row.read1), file(row.read2))}
        .set { samples }

    reference = file(params.reference_genome)
    reference_fna = file("ref/GCF_000149955.1_ASM14995v2_genomic.fna")
    annotation = file("ref/GCF_000149955.1_ASM14995v2_genomic.gff")
    inspectReads(samples, reference)
    fastqc(samples)

    preprocess_assembly(samples, reference)
    alignment_variant_calling(samples, reference)
    BLASTSearch(preprocess_assembly.out.sample_id, preprocess_assembly.out.chr0_contigs)

    liftoff(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold, reference_fna, annotation)
    dbscan(liftoff.out.sample_id, liftoff.out.genes_gff)
    blastp(dbscan.out.sample_id, dbscan.out.clustered_genes)
    antismash(liftoff.out.sample_id, liftoff.out.scaffold, liftoff.out.genes_gff)

    nucmerMummer(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold, reference_fna)
    repeatModeler(preprocess_assembly.out.sample_id, preprocess_assembly.out.scaffold)
    repeatMasker(repeatModeler.out.sample_id, preprocess_assembly.out.scaffold, repeatModeler.out.repeat_lib)
    // read1 = samples.map { it[1] }
    // read2 = samples.map { it[2] }
    // mcclintock(repeatMasker.out.sample_id, repeatMasker.out.masked_fasta, repeatModeler.out.repeat_lib, read1, read2, repeatMasker.out.masked_gff)
}
