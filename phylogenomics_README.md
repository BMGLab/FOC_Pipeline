# BUSCO-Based Phylogenomics Pipeline

This repository contains helper scripts for large-scale phylogenomic reconstruction of *Fusarium oxysporum* genomes using BUSCO single-copy orthologs. The workflow automates BUSCO parsing, ortholog extraction, multiple sequence alignment, concatenation, and phylogenetic tree inference with IQ-TREE and/or RAxML-NG. 

---

# Overview

The pipeline was developed for comparative phylogenomics of *Fusarium oxysporum* genomes, including both local isolates and publicly available assemblies. It supports:

* BUSCO-based single-copy ortholog extraction
* Automated MAFFT alignments
* Concatenated supermatrix generation
* Partition file creation
* Maximum-likelihood phylogenetic inference
* IQ-TREE3 and RAxML-NG support
* Filtering genomes using custom keep lists
* Reuse of precomputed alignments for faster reruns

---

# Scripts

## `fox_busco_phylogeny.py`

Main end-to-end phylogenomics workflow.

### Features

* Runs BUSCO on input genomes
* Detects existing BUSCO outputs
* Extracts single-copy BUSCO proteins
* Aligns orthologs using MAFFT
* Concatenates alignments into a supermatrix
* Generates partition files
* Runs IQ-TREE3 phylogenetic analyses

---

## `fox_busco_phylogeny2.py`

Modified version optimized for filtered datasets and reruns.

### Additional Features

* `--keep-samples` option for retaining selected genomes only
* Reuses existing alignments from `02_alignments`
* Skips repeated alignment generation
* Faster iterative phylogenetic reconstruction
* Includes optional `RAxML-NG` support

This version is particularly useful after Mash-based genome filtering or BUSCO completeness filtering. 

---

# Requirements

The following tools must be available in your environment:

* BUSCO
* MAFFT
* IQ-TREE3
* Python ≥ 3.9

Optional:

- RAxML-NG

Recommended conda environment:

```bash
conda activate fox-phylo
```

---

# Example command:

```bash
iqtree3 \
  -s 03_concatenated/fox_phylogenomics.concatenated.faa \
  -spp 03_concatenated/fox_phylogenomics.partitions.txt \
  -m LG+F+R \
  -tree-threads 16 \
  -iqtree-bootstrap 1000 \
  -seed 1337 \
  -safe \
  -pre 04_trees/iqtree/fox_phylogenomics


# Output Structure

```text
phylogenomics/
├── 01_single_copy_fastas/
├── 02_alignments/
├── 03_concatenated/
├── 04_trees/
└── metadata/
```

## Important Outputs

| File                 | Description                                      |
| -------------------- | ------------------------------------------------ |
| `*.concatenated.faa` | Concatenated BUSCO supermatrix                   |
| `*.partitions.txt`   | Partition coordinates for phylogenetic inference |
| `*.gene_stats.tsv`   | Per-gene alignment statistics                    |
| `*.sample_stats.tsv` | BUSCO gene counts per genome                     |
| `iqtree/`            | IQ-TREE outputs                                  |

---

# Phylogenetic Strategy

The workflow uses conserved single-copy BUSCO orthologs from the `hypocreales_odb10` lineage dataset to generate robust phylogenomic datasets suitable for fungal comparative genomics.

Default settings include:

* MAFFT `--auto`
* - IQ-TREE model: `LG+F+R`
* Ultrafast bootstrap: `1000`
* Random seed: `1337`

The pipeline automatically filters BUSCO genes according to minimum sample representation thresholds before concatenation. 

---

# Notes

* Existing BUSCO outputs are automatically reused unless `--force-busco` is specified.
* `fox_busco_phylogeny2.py` assumes alignments already exist in `02_alignments/`.
* The keep-list file should contain one sample ID per line.
* The workflow is optimized for large fungal genome collections.

