include { fastp } from '../modules/fastp.nf'
include { fastqc } from '../modules/fastqc.nf'
include { samTools } from '../modules/samtools.nf'

process megahitAssembly {
    tag "$id"
    publishDir "${params.output_dir}/assembly/${id}", mode: 'copy'

    input:
    tuple val(id), val(read1), val(read2)

    output:
    path "megahit_output/final.contigs.fa", emit: contigs
    path "megahit_output", emit: assembly_dir
    val id, emit: sample_id

    script:
    """
    source ${params.conda_shell}
    conda activate megahit
    megahit -1 '${read1}' -2 '${read2}' \
            -o 'megahit_output' \
            --presets meta-sensitive
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

process minimap2 {
    tag "$id"
    publishDir "${params.output_dir}/minimap2/${id}", mode: 'copy'

    input:
    val id
    path reference
    path scaffold

    output:
    tuple val(id), path("${id}_aligned.sam"), emit: sam

    script:
    """
    source $params.conda_shell
    conda activate minimap2
    minimap2 -ax sr ${reference} ${scaffold} > ${id}_aligned.sam
    conda deactivate
    """
}

process quast {
    tag "$id"
    publishDir "${params.output_dir}/assembly_qc/${id}", mode: 'copy'

    input:
    val id
    path contigs

    output:
    path "quast_output", emit: results

    script:
    """
    source $params.conda_shell
    export JAVA_HOME=\$HOME/miniconda3/envs/quast
    export JAVA_LD_LIBRARY_PATH=\${JAVA_LD_LIBRARY_PATH:-}
    conda activate quast
    quast.py ${contigs} \
             -o quast_output \
             --fast
    conda deactivate
    """
}

workflow preprocess_assembly {
    take:
    samples
    reference

    main:
    fastp(samples)
    fastqc(samples)
    megahit_out = megahitAssembly(fastp.out.trimmed_reads)
    quast(megahit_out.sample_id, megahit_out.contigs)
    ragtag_out = ragtagCorrect(megahit_out.sample_id, megahit_out.contigs, reference) | ragtagScaffold
    extractChr0Contigs(ragtag_out.sample_id, ragtag_out.scaffold_agp, ragtag_out.corrected_contigs)
    minimap2(ragtag_out.sample_id, ragtag_out.scaffold_fasta, reference)
    samTools(minimap2.out)

    emit:
    sample_id = ragtag_out.sample_id
    chr0_contigs = extractChr0Contigs.out.chr0_contigs
    scaffold = ragtag_out.scaffold_fasta
}