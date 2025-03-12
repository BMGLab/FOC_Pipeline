

process fastp {
    input:
    path(reads)
    
    output:
    path 'cleaned_reads.fastq.gz'
    
    script:
    """
    fastp -i $reads -o cleaned_reads.fastq.gz
    """
}

process megahit {
    input:
    path(reads)
    
    output:
    path 'assembly.fasta'
    
    script:
    """
    megahit -r cleaned_reads.fastq.gz -o megahit_output && mv megahit_output/final.contigs.fa assembly.fasta
    """
}

process ragtag_correction {
    input:
    path(assembly)
    path(reference)
    
    output:
    path 'ragtag.corrected.fasta'
    
    script:
    """
    ragtag.py correct -o ragtag_correction_output $reference $assembly
    mv ragtag_correction_output/ragtag.correct.fasta ragtag.corrected.fasta
    """
}

process ragtag_scaffold {
    input:
    path(corrected)
    path(reference)
    
    output:
    path ('ragtag.scaffolds.fasta'), emit: scaffolds_fasta
    path ('ragtag.scaffolds.agp'), emit: scaffolds_agp
    
    script:
    """
    ragtag.py scaffold -o ragtag_scaffold_output $reference $corrected
    mv ragtag_scaffold_output/ragtag.scaffolds.fasta ragtag.scaffolds.fasta
    mv ragtag_scaffold_output/ragtag.scaffolds.agp ragtag.scaffolds.agp
    """
}

process extract_chr0_contigs {
    input:
    path(agp)
    path(fasta)

    output:
    path 'chr0_contigs.fasta', emit: chr0_contigs

    script:
    """
    awk '\$5 == "U" {print \$6}' ${agp} | seqtk subseq ${fasta} - > chr0_contigs.fasta
    """
}

process nucmer_analysis {
    input:
    path(contigs)
    path(reference)
    
    output:
    path 'nucmer.delta'
    
    script:
    """
    nucmer --prefix=nucmer_output $reference $contigs
    """
}

process blastAnalysis {
    input:
    path(contigs)

    output:
    path 'blast_results.txt', emit: blast_results

    script:
    """
    blastn -query ${contigs} -db nt -out blast_results.txt -outfmt 6
    """
}

process repeatAnalysis {
    input:
    path(contigs)
    
    output:
    path 'repeatmasker.out'
    
    script:
    """
    RepeatMasker -species fungi ${contigs}
    """
}

workflow unplaced_contigs_analysis {
    reads_ch = Channel.fromPath(params.reads)
    reference_ch = Channel.fromPath(params.reference)

    fastp(reads_ch) | megahit 
    ragtag_correction(megahit.out, reference_ch)
    ragtag_scaffold(ragtag_correction.out, reference_ch)
    ragtag_scaffold.out
    extract_chr0_contigs(ragtag_scaffold.out.scaffolds_agp, ragtag_scaffold.out.scaffolds_fasta)
    nucmer_analysis(extract_chr0_contigs.out.chr0_contigs, reference_ch)
    blastAnalysis(extract_chr0_contigs.out.chr0_contigs)
    repeatAnalysis(extract_chr0_contigs.out.chr0_contigs)
}
