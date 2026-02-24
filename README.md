# Whole-genome assemblies and comparative genomics of Fusarium oxysporum f. sp. capsici Pipeline

Nextflow DSL2 pipeline for assembly/scaffolding, coverage and variant analysis, repeat profiling, and downstream protein/functional analyses.

## Entry points

- `main.nf`
- `workflows/preprocess_assembly.nf`
- `workflows/alignment_variant_calling.nf`

## Workflow summary

1. Read QC and trimming
- `FastQC`, `fastp`

2. Assembly and scaffolding
- `megahit`, `ragtag`, `QUAST`, `LAST`
- Optional Chr0 append flow via BLAST (`BLASTn` + `appendContigs`)

3. Coverage estimation
- Read-to-reference mapping with `bwa mem`
- Duplicate handling and depth/breadth summary using `samtools`
- Per-sample metrics plus combined coverage summary

4. Variant calling and annotation
- `Bowtie2` alignment
- `bcftools mpileup/call/filter` with configurable MQ/BQ/ploidy thresholds
- Callable-region workflow: `mosdepth` + `samtools depth`
- Consensus callable BED generation across samples
- Restriction of variants to callable regions
- `SnpEff` annotation

5. Group-specific variant analysis (optional)
- Strict group-specific variant extraction from callable-filtered VCFs
- Optional downstream annotation branch:
  - Group VCF annotation
  - Candidate gene/protein extraction
  - `PHI-base` BLASTp
  - `EggNOG` annotation on Fol4287 proteins

6. Repeat/TE and protein analyses
- `RepeatModeler`, `RepeatMasker`
- TE summary tables per sample
- Optional `McClintock`
- `Liftoff` gene transfer + protein extraction
- `SIX` BLASTp screening
- `dbCAN` (HMM mode)
- `TargetP`, `SignalP`, `WoLFPSort`

## Inputs

### Samplesheet
Default: `samplesheet.csv`

```csv
sample_id,read1,read2
34,/path/sample34_R1.fastq.gz,/path/sample34_R2.fastq.gz
35,/path/sample35_R1.fastq.gz,/path/sample35_R2.fastq.gz
```

### Group map
Used when `--avc_run_group_specific true`.
Default: `samplesheets/group_map.tsv`

```tsv
# sample_id	group
34	GroupB
35	GroupB
41	GroupA
42	GroupA
```

### External resources
Configure in `nextflow.config` or via CLI params:

- BLAST DBs and query FASTA (including SIX queries)
- dbCAN database
- PHI-base BLAST database
- EggNOG data directory
- SnpEff DB/JAR paths
- Tool paths and conda shell path

## Run

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/sercanozturk/miniconda3/envs/nf-env
nextflow run main.nf --samplesheet_path samplesheet.csv -with-report
```

## Key parameters

- `--output_dir`
- `--samplesheet_path`
- `--reference_genome`
- `--enable_chr0_blast_append`
- `--run_mcclintock`
- `--coverage_threads`
- `--six_queries_fasta`
- `--avc_group_map`
- `--avc_run_group_specific`
- `--avc_run_group_functional_annotation`
- `--avc_min_mapq`
- `--avc_min_baseq`
- `--avc_min_callable_depth`
- `--avc_consensus_callable_fraction`
- `--avc_ploidy`

## Outputs (high-level)

Default base: `results/`

- `coverage/`
- `assembly/`, `ragtag/`, `chr0_contigs/`, `appended_contigs/`
- `alignment/`, `samtools/`, `bcftools/`, `bcftools_callable/`, `snpeff/`, `callable/`
- `group_specific_variants/`, `eggnog/`
- `repeatmodeler/`, `repeatmasker/`, `te_summary/`, `mcclintock/`
- `liftoff/`, `dbcan/`, `six_blastp/`, `targetp/`, `signalp/`, `wolfpsort/`

## Notes

- Optional branches are controlled by config flags.
- Validate all database/tool paths before production runs.
