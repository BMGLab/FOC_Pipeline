# Comparative Genomic Analysis of Fusarium oxysporum f. sp. lycopercisi Pipeline

This repository contains the Nextflow DSL2 workflow used for comparative genomic analyses of local *Fusarium oxysporum* isolates together with publicly available genomes. The pipeline integrates genome assembly, scaffolding, repeat masking, gene transfer, variant calling, secretome prediction, effector analysis, and functional annotation.


## What this pipeline does

Given paired-end Illumina reads and a reference genome, the workflow performs:

- Read quality control and trimming
- De novo genome assembly and reference-guided scaffolding
- Whole-genome alignment and assembly quality evaluation
- Repeat and transposable element discovery
- Reference-guided gene transfer and protein extraction
- Variant calling and functional impact prediction
- Secretome and effector prediction
- CAZyme and functional annotation
- SIX gene screening

The workflow was primarily developed for analyses of local *Fusarium oxysporum f. sp. lycopersici* (Fol) isolates using the Fol4287 reference genome.

---

# Main workflow structure

Main entrypoint:

```text
main.nf

Sub-workflows:
- `workflows/preprocess_assembly.nf`
- `workflows/alignment_variant_calling.nf`

Additional modules and utility scripts are organized under:
modules/
utilities/
```

# Workflow summary
1. Read QC and trimming:Raw Illumina paired-end reads are filtered and quality-controlled before downstream analyses.

### Tools

-FastQC
-fastp

### Outputs

```text
fastqc/
fastp/
```
---

2. Genome assembly and scaffolding:Reads are assembled de novo and scaffolded against the Fol4287 reference genome.

### Tools

-MEGAHIT
-RagTag
-QUAST
-LAST
-nucmer
-mummerplot

### Workflow

```text
fastp → MEGAHIT → RagTag → QUAST → LAST/nucmer
```

### Outputs

```text
assembly/
ragtag/
chr0_contigs/
quast/
synteny/
```

---

3. Repeat and transposable element analysis:Repeat regions and transposable elements are identified and summarized across genomes.

### Tools

-RepeatModeler
-RepeatMasker
-McClintock (optional)

### Workflow

```text
RepeatModeler → RepeatMasker → TE analysis
```

### Outputs

```text
repeatmodeler/
repeatmasker/
mcclintock/
te_summary/
```

---

4. Gene transfer and annotation:Reference-guided gene transfer and protein extraction are performed using the Fol4287 annotation.

### Tools 

-Liftoff
-AUGUSTUS
-AGAT

### Outputs

```text
liftoff/
augustus/
proteins/
```

5. Variant calling and functional impact analysis:Reference-based variant analyses are performed by mapping reads to the Fol4287 reference genome.

### Tools

- Bowtie2
- samtools
- bcftools
- mosdepth
- SnpEff

### Workflow

```text
Bowtie2 → samtools → bcftools → callable filtering → SnpEff
```

### Outputs

```text
alignment/
samtools/
bcftools/
coverage/
callable/
snpeff/
```

---

6. Functional annotation and secretome analysis:Predicted proteins are analyzed to identify secreted proteins, candidate effectors, and carbohydrate-active enzymes.


### Tools

- SignalP
- TargetP
- WoLFPSort
- EffectorP
- dbCAN
- antiSMASH
- BLASTp vs MEROPS
- eggNOG-mapper
- DeepTMHMM (optional)

### Workflow

```text
Protein FASTA → SignalP/TargetP/WoLFPSort → secretome prediction
Protein FASTA → EffectorP → effector prediction
Protein FASTA → dbCAN → CAZyme annotation
Protein FASTA → antiSMASH → secondary metabolite prediction
Protein FASTA → BLASTp vs MEROPS → protease annotation
```

### Outputs

```text
signalp/
targetp/
wolfpsort/
effectorp/
dbcan/
antismash/
blastp/
eggnog/
deeptmhmm/
```

7. SIX gene analysis:SIX homologs are identified using BLAST-based homology searches against genome assemblies.

### Tools

- BLAST+
- tblastn

### Workflow

```text
SIX queries → tblastn → filtering → presence/absence matrix
```


# Repository layout

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
- The workflow was developed for comparative genomics analyses of Fusarium oxysporum isolates.
- Fol4287 was used as the primary reference genome.
- BUSCO-based phylogenomics analyses are documented separately.
- Multiple conda environments are used internally by different workflow processes.
- Optional branches are controlled through parameters in nextflow.config.
- `workflows/functional_annotation.nf` is present but not called by `main.nf`.
- Slack notification is sent on completion if `SLACK_WEBHOOK` is set.
