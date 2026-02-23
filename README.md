# Whole-genome assemblies and comparative genomics of Fusarium oxysporum f. sp. capsici Pipeline

Nextflow DSL2 pipeline for assembly processing, variant discovery, TE analysis, gene/protein functional annotation, and group-specific variant interpretation.

## Main entrypoints

- `main.nf`
- `workflows/preprocess_assembly.nf`
- `workflows/alignment_variant_calling.nf`

## Pipeline components

- Read QC and trimming: `FastQC`, `fastp`
- Assembly/scaffolding: `megahit`, `ragtag`
- Coverage estimation: `bwa mem` + `samtools` depth metrics
- Variant calling: `Bowtie2`, `bcftools`, `SnpEff`
- Callable genome workflow: `mosdepth` + `samtools depth` + consensus callable BED
- Group-specific variants: strict per-group variant extraction
- TE analysis: `RepeatModeler`, `RepeatMasker` (optional `McClintock`)
- Gene/protein analyses: `Liftoff`, `antiSMASH`, `BLASTp` (MEROPS), `dbCAN`, `SIX BLASTp`, `TargetP`, `SignalP`, `WoLFPSort`, optional `DeepTMHMM`
- Functional annotation for group-specific candidates: EggNOG + PHI-base BLASTp

## Required inputs

### 1) Samplesheet

CSV with header:

```csv
sample_id,read1,read2
34,/path/sample34_R1.fastq.gz,/path/sample34_R2.fastq.gz
35,/path/sample35_R1.fastq.gz,/path/sample35_R2.fastq.gz
```

Default: `samplesheet.csv`

### 2) Group map

TSV format:

```tsv
# sample_id	group
34	GroupB
35	GroupB
41	GroupA
42	GroupA
```

Default: `samplesheets/group_map.tsv`

### 3) External resources

Configure paths in `nextflow.config` or via CLI:

- BLAST databases (nt/MEROPS/PHI-base)
- SIX query FASTA
- dbCAN database
- EggNOG data directory
- SnpEff DB/JAR
- Tool paths (Conda env/tool binaries)

## Running

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/sercanozturk/miniconda3/envs/nf-env
nextflow run main.nf --samplesheet_path samplesheet.csv -with-report
```

## Key parameters

- `--output_dir`
- `--samplesheet_path`
- `--avc_group_map`
- `--avc_run_group_specific`
- `--avc_run_group_functional_annotation`
- `--enable_chr0_blast_append`
- `--run_mcclintock`
- `--skip_deeptmhmm`

## Outputs (high level)

Under `results/` (default):

- `coverage/`
- `callable/`
- `bcftools_callable/`
- `group_specific_variants/`
- `eggnog/`
- `six_blastp/`
- `dbcan/`
- `repeatmodeler/`, `repeatmasker/`
- `liftoff/`, `antismash/`, `blastp/`, `snpeff/`

## Notes

- Optional analysis modules can be enabled/disabled with config flags.
- Validate group map and database paths before production runs.
