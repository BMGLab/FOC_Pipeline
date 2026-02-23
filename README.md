# FoxGenome_wd (`dogay` branch)

This branch contains the manuscript-aligned implementation of the *Fusarium oxysporum* pipeline.

## Scope of this branch

The `dogay` branch is configured to match the manuscript workflow as closely as possible, including:

- Read mapping for coverage estimation with `bwa mem`, duplicate handling, and depth-based summary.
- RagTag assembly/scaffolding as core assembly step.
- Variant calling with `Bowtie2 + bcftools` using explicit ploidy and quality thresholds.
- Callable genome workflow (`mosdepth` + `samtools depth`) with consensus callable genome (K-of-N, default 80%).
- Strict group-specific variant extraction from callable-restricted variants.
- Functional annotation of group-specific variants (SnpEff, EggNOG on Fol4287 proteins, PHI-base BLASTp).
- SIX effector screening using dedicated BLASTp process and manuscript thresholds.
- CAZyme identification in HMM-only mode with `domE <= 1e-15` and protein-level summaries.

Optional/non-manuscript analyses are retained but disabled by default where applicable.

## Main entrypoints

- `main.nf`
- `workflows/preprocess_assembly.nf`
- `workflows/alignment_variant_calling.nf`

## Key manuscript-alignment defaults

From `nextflow.config` (can be overridden at runtime):

- `params.enable_chr0_blast_append = false`
- `params.run_mcclintock = false`
- `params.avc_min_mapq = 30`
- `params.avc_min_baseq = 20`
- `params.avc_min_callable_depth = 20`
- `params.avc_consensus_callable_fraction = 0.80`
- `params.avc_ploidy = 2`
- `params.avc_run_group_specific = true`
- `params.avc_run_group_functional_annotation = true`
- `params.dbcan_hmm_dome = "1e-15"`

## Required inputs

### 1) Samplesheet

CSV with header:

```csv
sample_id,read1,read2
34,/path/sample34_R1.fastq.gz,/path/sample34_R2.fastq.gz
35,/path/sample35_R1.fastq.gz,/path/sample35_R2.fastq.gz
```

Default: `samplesheet.csv`

### 2) Group map for strict group-specific variants

TSV with header or comment lines allowed:

```tsv
# sample_id	group
34	GroupB
35	GroupB
41	GroupA
42	GroupA
```

Default: `samplesheets/group_map.tsv`

### 3) External databases/resources

Set valid paths in `nextflow.config` or via CLI params:

- BLAST nt/MEROPS and SIX query fasta
- dbCAN database directory
- PHI-base BLAST database (`params.phibase_db`)
- EggNOG data directory (`params.eggnog_data_dir`)
- SnpEff DB/JAR paths

## Running (dogay branch)

Use your validated environment:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/sercanozturk/miniconda3/envs/nf-env
nextflow run main.nf \
  --samplesheet_path samplesheet.csv \
  -with-report
```

## Outputs (high level)

Under `results/` (default):

- `coverage/` read-based coverage metrics
- `callable/` per-sample and consensus callable BEDs
- `bcftools_callable/` variants restricted to callable regions
- `group_specific_variants/` strict group-specific VCFs + summaries
- `group_specific_variants/phi_base/` PHI-base hits
- `eggnog/` Fol4287 EggNOG annotations
- `six_blastp/` SIX screening results
- `dbcan/` HMM-filtered and protein-level CAZyme summaries

## Notes

- This branch keeps several legacy/extended processes for flexibility, but manuscript-aligned path is enabled by default.
- Ensure sample groups and database paths are correct before production runs.
