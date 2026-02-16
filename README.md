# FoxGenome_wd

Nextflow DSL2 workflow for fungal genome assembly refinement, variant calling, repeat/transposon analysis, and functional annotation.

## What this pipeline does

Given paired-end reads and a reference genome, the pipeline runs:

1. Read QC and trimming (`FastQC`, `fastp`)
2. Assembly and scaffolding (`megahit`, `ragtag`)
3. Contig handling around `Chr0_RagTag` and scaffold merging
4. Alignment and variant calling (`bowtie2`, `samtools`, `bcftools`, `SnpEff`)
5. Assembly comparison and quality checks (`nucmer`, `mummerplot`, `QUAST`, `LAST`)
6. Repeat discovery/masking and TE analysis (`RepeatModeler`, `RepeatMasker`, `McClintock`)
7. Gene prediction/transfer (`AUGUSTUS`, `Liftoff`)
8. Functional analyses (`antiSMASH`, `BLASTp` vs MEROPS, `dbCAN`, `TargetP`, `SignalP`, `WoLFPSort`, optional `DeepTMHMM`)

Main entrypoint: `main.nf`  
Sub-workflows:
- `workflows/preprocess_assembly.nf`
- `workflows/alignment_variant_calling.nf`

## Repository layout

```text
FoxGenome_wd/
├── main.nf
├── nextflow.config
├── workflows/
├── modules/
├── utilities/
├── ref/
├── sample_reads/
├── samplesheet.csv
└── run.sh
```

## Requirements

- Nextflow (tested with `NXF_VER=24.10.5` in `run.sh`)
- Conda installation with tool-specific environments used by the workflow
- Local databases and external resources referenced in `nextflow.config`

Important: this project expects existing local paths for some resources (example: BLAST db, dbCAN db, McClintock, TargetP). Update these in `nextflow.config` for your system:

- `params.blast_db_dir`
- `params.dbcan_db_dir`
- `params.mcclintock`
- `params.targetp`
- `params.merops_db_dir`
- `params.conda_shell`

## Input format

Samplesheet CSV (header required):

```csv
sample_id,read1,read2
sampleA,/abs/path/sampleA_R1.fastq.gz,/abs/path/sampleA_R2.fastq.gz
sampleB,/abs/path/sampleB_R1.fastq.gz,/abs/path/sampleB_R2.fastq.gz
```

Default samplesheet path is `samplesheet.csv` (override with `--samplesheet_path`).

## Running

Example run:

```bash
export NXF_VER=24.10.5
nextflow -log logs/.nextflow.log run main.nf \
  --samplesheet_path samplesheet.csv \
  --skip_deeptmhmm true
```

The provided `run.sh` contains a similar command and optional Slack webhook setup.

## Key parameters

- `--output_dir` (default: `results`)
- `--reference_genome` (default: `ref/GCF_000149955.1_ASM14995v2_genomic.fna.gz`)
- `--samplesheet_path` (default: `samplesheet.csv`)
- `--skip_deeptmhmm` (default: `false`)

See `nextflow.config` for additional tool/database settings.

## Outputs

Results are published under `${output_dir}` (default `results`) by analysis type, including:

- `fastp/`, `fastqc/`
- `assembly/`, `ragtag/`, `chr0_contigs/`
- `alignment/`, `samtools/`, `bcftools/`, `snpeff/`
- `repeatmodeler/`, `repeatmasker/`, `mcclintock/`
- `liftoff/`, `augustus/`, `antismash/`, `blastp/`, `dbcan/`
- `targetp/`, `signalp/`, `wolfpsort/`, `deeptmhmm/` (if enabled)

## Notes

- `workflows/functional_annotation.nf` is present but not called by `main.nf`.
- The pipeline uses multiple conda environments inside process scripts.
- Slack notification is sent on completion if `SLACK_WEBHOOK` is set.
