#!/usr/bin/env python3
"""Fox genome phylogenomics helper.

This script automates the pipeline:
  genome FASTA -> BUSCO -> single-copy extraction -> alignment/concatenation -> IQ-TREE/RAxML.

Default paths assume the current FoxGenome_wd checkout layout, but most CLI
options can be overridden to operate on any dataset.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import logging
import re
import shutil
import subprocess
import sys
from collections import defaultdict, OrderedDict
from dataclasses import dataclass
from math import ceil
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

LOGGER = logging.getLogger("fox_busco_phylogeny")


@dataclass
class BuscoRun:
    sample_id: str
    run_dir: Path

    @property
    def single_copy_dir(self) -> Path:
        return self.run_dir / "busco_sequences" / "single_copy_busco_sequences"


def guess_project_root() -> Path:
    return Path(__file__).resolve().parent.parent


PROJECT_ROOT = guess_project_root()
DEFAULT_PUBLIC_BUSCO = PROJECT_ROOT / "beyza_article_analysis/ncbi_genome_comparison/busco_results"
DEFAULT_FOL_BUSCO = PROJECT_ROOT / "beyza_article_analysis/ncbi_genome_comparison/fol_genomes"
DEFAULT_OUT_DIR = PROJECT_ROOT / "beyza_article_analysis/ncbi_genome_comparison/phylogenomics"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run BUSCO-driven phylogenomics on Fox genomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--genome-dir", type=Path, help="Directory with genome FASTA files", default=None)
    parser.add_argument(
        "--genome-patterns",
        nargs="+",
        default=("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fa.gz", "*.fasta.gz"),
        help="Glob patterns (relative to --genome-dir) to pick genome files",
    )
    parser.add_argument("--genome-list", type=Path, default=None, help="Optional text file with one genome FASTA path per line")
    parser.add_argument("--busco-bin", default="busco", help="BUSCO executable")
    parser.add_argument("--busco-lineage", default="hypocreales_odb10", help="BUSCO lineage dataset")
    parser.add_argument("--busco-mode", choices=("genome", "proteins", "transcriptome"), default="genome")
    parser.add_argument("--busco-threads", type=int, default=16, help="Threads per BUSCO job")
    parser.add_argument("--busco-workers", type=int, default=2, help="Parallel BUSCO jobs to run at once")
    parser.add_argument("--busco-out", type=Path, default=DEFAULT_FOL_BUSCO / "busco_custom_runs", help="Directory to store new BUSCO outputs")
    parser.add_argument("--busco-extra", nargs="*", default=(), help="Extra arguments appended to BUSCO")
    parser.add_argument("--skip-busco", action="store_true", help="Skip BUSCO runs even if genomes are supplied")
    parser.add_argument("--force-busco", action="store_true", help="Re-run BUSCO even if outputs already exist")
    parser.add_argument(
        "--busco-roots",
        nargs="+",
        default=(DEFAULT_FOL_BUSCO, DEFAULT_PUBLIC_BUSCO),
        help="Directories that already contain BUSCO run folders",
    )

    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Pipeline working/output directory")
    parser.add_argument("--prefix", default="fox_phylogenomics", help="Prefix for concatenated alignment + tree outputs")
    parser.add_argument("--min-sample-fraction", type=float, default=0.6, help="Minimum fraction of samples per gene")
    parser.add_argument("--min-samples", type=int, default=4, help="Absolute minimum samples per gene (overrides fraction if larger)")

    parser.add_argument("--aligner", choices=("mafft", "none"), default="mafft", help="Alignment engine for BUSCO sequences")
    parser.add_argument("--mafft-bin", default="mafft", help="MAFFT executable path")
    parser.add_argument("--mafft-threads", type=int, default=4, help="Threads per MAFFT alignment")
    parser.add_argument("--mafft-extra", nargs="*", default=("--auto",), help="Extra MAFFT args")

    parser.add_argument(
        "--tree-programs",
        nargs="+",
        choices=("iqtree3", "raxml-ng", "none"),
        default=("iqtree3",),
        help="Tree builders to run after concatenation",
    )
    parser.add_argument("--iqtree-bin", default="iqtree3", help="IQ-TREE executable")
    parser.add_argument("--iqtree-model", default="LG+F+I+G4", help="IQ-TREE substitution model or MFP directive")
    parser.add_argument("--iqtree-bootstrap", type=int, default=1000, help="Ultrafast bootstrap replicates (-bb)")
    parser.add_argument("--iqtree-extra", nargs="*", default=(), help="Extra IQ-TREE args")

    parser.add_argument("--raxml-bin", default="raxml-ng", help="RAxML-NG executable")
    parser.add_argument("--raxml-model", default="LG+G4", help="RAxML-NG model (or partition file)")
    parser.add_argument("--raxml-bootstrap", type=int, default=200, help="RAxML-NG bootstrap trees")
    parser.add_argument("--raxml-extra", nargs="*", default=(), help="Extra RAxML-NG args")

    parser.add_argument("--tree-threads", type=int, default=32, help="Threads for tree inference")
    parser.add_argument("--seed", type=int, default=1337, help="Random seed for reproducibility")

    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing heavy steps")
    parser.add_argument("--log-level", default="INFO", help="Logging level (DEBUG, INFO, WARNING...")
    parser.add_argument(
    "--keep-samples",
    type=Path,
    default=None,
    help="Optional text file with sample IDs to keep (one per line)"
)

    return parser.parse_args(argv)


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(asctime)s] %(levelname)s - %(message)s",
    )


def ensure_executable(exe: str) -> None:
    if shutil.which(exe) is None:
        raise FileNotFoundError(f"Executable '{exe}' not found in PATH")


def collect_genomes(genome_dir: Path | None, patterns: Sequence[str], genome_list: Path | None) -> List[Path]:
    genomes: Dict[Path, None] = {}
    if genome_dir:
        for pattern in patterns:
            for path in genome_dir.glob(pattern):
                if path.is_file():
                    genomes[path.resolve()] = None
    if genome_list and genome_list.exists():
        with genome_list.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                path = Path(line).expanduser().resolve()
                if path.is_file():
                    genomes[path] = None
                else:
                    LOGGER.warning("Genome path %s from %s does not exist", path, genome_list)
    return sorted(genomes.keys())


def infer_sample_id_from_genome(genome_path: Path) -> str:
    name = genome_path.name
    name = re.sub(r"\.(fasta|fa|fna)(\.gz)?$", "", name, flags=re.IGNORECASE)
    return name


def busco_already_done(out_dir: Path, sample_id: str, lineage: str) -> bool:
    run_dir = out_dir / sample_id / f"run_{lineage}"
    single_dir = run_dir / "busco_sequences" / "single_copy_busco_sequences"
    return single_dir.exists()


def run_busco_job(
    genome_path: Path,
    args: argparse.Namespace,
) -> Tuple[str, int]:
    sample_id = infer_sample_id_from_genome(genome_path)
    out_dir = args.busco_out / sample_id
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / f"busco_{sample_id}.log"
    cmd = [
        args.busco_bin,
        "-i",
        str(genome_path),
        "-l",
        args.busco_lineage,
        "-o",
        sample_id,
        "-m",
        args.busco_mode,
        "-c",
        str(args.busco_threads),
        "--out_path",
        str(args.busco_out),
        "--offline",
        "--download_path",
        str(args.busco_out / "busco_downloads"),
    ]
    if args.force_busco:
        cmd.append("--force")
    if args.busco_extra:
        cmd.extend(args.busco_extra)
    LOGGER.info("BUSCO: %s", " ".join(cmd))
    if args.dry_run:
        return sample_id, 0
    ensure_executable(args.busco_bin)
    args.busco_out.mkdir(parents=True, exist_ok=True)
    args.busco_out.joinpath("busco_downloads").mkdir(exist_ok=True)
    with log_path.open("w") as log_handle:
        proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT)
    return sample_id, proc.returncode


def run_busco_stage(genomes: Sequence[Path], args: argparse.Namespace) -> None:
    if not genomes:
        LOGGER.info("No genomes provided for BUSCO stage")
        return
    skipped = 0
    jobs: List[Path] = []
    for genome in genomes:
        sample_id = infer_sample_id_from_genome(genome)
        if not args.force_busco and busco_already_done(args.busco_out, sample_id, args.busco_lineage):
            LOGGER.info("Skipping BUSCO for %s (already exists)", sample_id)
            skipped += 1
            continue
        jobs.append(genome)
    if not jobs:
        LOGGER.info("All genomes already processed by BUSCO")
        return
    LOGGER.info("Running BUSCO on %d genomes (%d skipped)", len(jobs), skipped)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.busco_workers) as executor:
        futures = [executor.submit(run_busco_job, genome, args) for genome in jobs]
        for future in concurrent.futures.as_completed(futures):
            sample_id, code = future.result()
            if code != 0:
                LOGGER.error("BUSCO failed for %s (exit %s)", sample_id, code)
                raise RuntimeError(f"BUSCO failed for {sample_id}")
    LOGGER.info("BUSCO stage completed")


def infer_sample_id_from_run(parent_name: str) -> str:
    return re.sub(r"_busco$", "", parent_name)


def discover_busco_runs(roots: Iterable[Path], lineage: str) -> Dict[str, BuscoRun]:
    runs: Dict[str, BuscoRun] = OrderedDict()
    pattern = f"run_{lineage}"
    for root in roots:
        root = Path(root)
        if not root.exists():
            LOGGER.warning("BUSCO root %s does not exist; skipping", root)
            continue
        for run_dir in root.rglob(pattern):
            if not run_dir.is_dir():
                continue
            sample_id = infer_sample_id_from_run(run_dir.parent.name)
            single_dir = run_dir / "busco_sequences" / "single_copy_busco_sequences"
            if not single_dir.exists():
                continue
            if sample_id in runs:
                LOGGER.debug("Sample %s already recorded; keeping first occurrence", sample_id)
                continue
            runs[sample_id] = BuscoRun(sample_id=sample_id, run_dir=run_dir)
    LOGGER.info("Discovered %d BUSCO runs", len(runs))
    return runs


def build_gene_index(busco_runs: Dict[str, BuscoRun]) -> Dict[str, Dict[str, Path]]:
    gene_index: Dict[str, Dict[str, Path]] = defaultdict(dict)
    for sample_id, run in busco_runs.items():
        single_dir = run.single_copy_dir
        for fasta in single_dir.glob("*.faa"):
            gene_id = fasta.stem
            gene_index[gene_id][sample_id] = fasta
    LOGGER.info("Indexed %d BUSCO genes", len(gene_index))
    return gene_index


def write_fasta_record(handle, header: str, seq: str, width: int = 60) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def read_single_fasta_sequence(path: Path) -> str:
    seq: List[str] = []
    with path.open() as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq)


def emit_gene_fastas(
    gene_index: Dict[str, Dict[str, Path]],
    samples: Sequence[str],
    min_required: int,
    dest_dir: Path,
) -> List[Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for gene_id in sorted(gene_index):
        sample_map = gene_index[gene_id]
        if len(sample_map) < min_required:
            continue
        out_path = dest_dir / f"{gene_id}.faa"
        with out_path.open("w") as handle:
            for sample_id in sorted(sample_map):
                seq = read_single_fasta_sequence(sample_map[sample_id])
                write_fasta_record(handle, f"{sample_id}|{gene_id}", seq)
        written.append(out_path)
    LOGGER.info("Wrote %d per-gene FASTA files", len(written))
    return written


def run_mafft(input_fasta: Path, output_fasta: Path, args: argparse.Namespace) -> None:
    cmd = [args.mafft_bin, "--thread", str(args.mafft_threads)] + list(args.mafft_extra) + [str(input_fasta)]
    LOGGER.debug("MAFFT: %s", " ".join(cmd))
    ensure_executable(args.mafft_bin)
    with output_fasta.open("w") as out_handle:
        proc = subprocess.run(cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        LOGGER.error("MAFFT failed on %s: %s", input_fasta.name, proc.stderr)
        raise RuntimeError(f"MAFFT failed for {input_fasta}")


def align_genes(per_gene_fastas: Sequence[Path], align_dir: Path, args: argparse.Namespace) -> List[Path]:
    align_dir.mkdir(parents=True, exist_ok=True)
    aligned: List[Path] = []
    for fasta in per_gene_fastas:
        out_path = align_dir / (fasta.stem + ".aln.faa")
        if args.aligner == "none":
            shutil.copy2(fasta, out_path)
        else:
            if args.dry_run:
                LOGGER.info("DRY-RUN MAFFT %s -> %s", fasta.name, out_path.name)
            else:
                run_mafft(fasta, out_path, args)
        aligned.append(out_path)
    LOGGER.info("Generated %d aligned genes", len(aligned))
    return aligned


def concatenate_alignments(
    aligned_fastas: Sequence[Path],
    samples: Sequence[str],
    out_dir: Path,
    prefix: str,
) -> Tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    concat_path = out_dir / f"{prefix}.concatenated.faa"
    partition_path = out_dir / f"{prefix}.partitions.txt"
    gene_stats_path = out_dir / f"{prefix}.gene_stats.tsv"

    sample_sequences: Dict[str, List[str]] = {sample: [] for sample in samples}
    partitions: List[Tuple[str, int, int]] = []
    with gene_stats_path.open("w") as stats_handle:
        stats_handle.write("gene_id\tlength\tn_samples\tsamples\n")
        cursor = 1
        for aln in sorted(aligned_fastas):
            records = parse_fasta(aln)
            if not records:
                continue
            lengths = {len(seq) for seq in records.values()}
            if len(lengths) != 1:
                raise ValueError(f"Alignment {aln} has sequences with inconsistent lengths")
            aln_len = lengths.pop()
            gene_id = aln.stem.replace(".aln", "")
            for sample in samples:
                seq = records.get(sample)
                if seq is None:
                    seq = "-" * aln_len
                sample_sequences[sample].append(seq)
            start = cursor
            end = cursor + aln_len - 1
            partitions.append((gene_id, start, end))
            stats_handle.write(
                f"{gene_id}\t{aln_len}\t{len(records)}\t{','.join(sorted(records))}\n"
            )
            cursor = end + 1
    with concat_path.open("w") as concat_handle:
        for sample in samples:
            seq = "".join(sample_sequences[sample])
            write_fasta_record(concat_handle, sample, seq)
    with partition_path.open("w") as part_handle:
        for gene_id, start, end in partitions:
            part_handle.write(f"LG, {gene_id} = {start}-{end}\n")
    return concat_path, partition_path, gene_stats_path


def parse_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, str] = OrderedDict()
    current_id = None
    seq_lines: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(seq_lines)
                header = line[1:]
                current_id = header.split("|")[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if current_id is not None:
            records[current_id] = "".join(seq_lines)
    return records


def run_iqtree(concat_fasta: Path, partition_file: Path, out_dir: Path, args: argparse.Namespace) -> None:
    ensure_executable(args.iqtree_bin)
    cmd = [
        args.iqtree_bin,
        "-s",
        str(concat_fasta),
        "-spp",
        str(partition_file),
        "-m",
        args.iqtree_model,
        "-nt",
        str(args.tree_threads),
        "-bb",
        str(args.iqtree_bootstrap),
        "-seed",
        str(args.seed),
        "-pre",
        str(out_dir / "iqtree" / args.prefix),
    ]
    if args.iqtree_extra:
        cmd.extend(args.iqtree_extra)
    out_tree_dir = out_dir / "iqtree"
    out_tree_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Running IQ-TREE: %s", " ".join(cmd))
    if args.dry_run:
        return
    subprocess.run(cmd, check=True)


def run_raxml(concat_fasta: Path, partition_file: Path, out_dir: Path, args: argparse.Namespace) -> None:
    ensure_executable(args.raxml_bin)
    model_arg = str(partition_file) if args.raxml_model == "partition" else args.raxml_model
    cmd = [
        args.raxml_bin,
        "--all",
        "--msa",
        str(concat_fasta),
        "--model",
        model_arg,
        "--threads",
        str(args.tree_threads),
        "--seed",
        str(args.seed),
        "--bs-trees",
        str(args.raxml_bootstrap),
        "--prefix",
        str(out_dir / "raxml" / args.prefix),
    ]
    if args.raxml_extra:
        cmd.extend(args.raxml_extra)
    out_tree_dir = out_dir / "raxml"
    out_tree_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Running RAxML-NG: %s", " ".join(cmd))
    if args.dry_run:
        return
    subprocess.run(cmd, check=True)


def summarize_samples(gene_index: Dict[str, Dict[str, Path]], samples: Sequence[str], out_dir: Path, prefix: str) -> Path:
    stats_path = out_dir / f"{prefix}.sample_stats.tsv"
    with stats_path.open("w") as handle:
        handle.write("sample_id\tsingle_copy_genes\n")
        for sample in samples:
            count = sum(1 for gene in gene_index.values() if sample in gene)
            handle.write(f"{sample}\t{count}\n")
    return stats_path


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    configure_logging(args.log_level)

    args.out_dir = args.out_dir.expanduser().resolve()
    args.busco_out = args.busco_out.expanduser().resolve()
    args.busco_roots = [Path(p).expanduser().resolve() for p in args.busco_roots]

    genomes = collect_genomes(args.genome_dir, args.genome_patterns, args.genome_list)
    if genomes and args.skip_busco:
        LOGGER.warning("--skip-busco set but genomes were supplied; genomes will be ignored for BUSCO stage")
    if genomes and not args.skip_busco:
        LOGGER.info("Starting BUSCO on %d genomes", len(genomes))
        run_busco_stage(genomes, args)
    elif genomes:
        LOGGER.info("BUSCO stage skipped (per user request)")

    normalized_roots: List[Path] = []
    if args.busco_out.exists() and args.busco_out not in args.busco_roots:
        normalized_roots.append(args.busco_out)
    for root in args.busco_roots:
        if root not in normalized_roots:
            normalized_roots.append(root)
    args.busco_roots = normalized_roots

    busco_runs = discover_busco_runs(args.busco_roots, args.busco_lineage)
    if not busco_runs:
        LOGGER.error("No BUSCO runs found for lineage %s", args.busco_lineage)
        return 1

# -------- KEEP FILTER --------
    if args.keep_samples:
        keep_ids = {
            line.strip().replace(">", "")
            for line in args.keep_samples.open()
            if line.strip()
        }

        original_count = len(busco_runs)

        busco_runs = {
            sid: run
            for sid, run in busco_runs.items()
            if sid in keep_ids or any(sid.startswith(k) for k in keep_ids)
        }

        LOGGER.info(
            "Filtered samples: kept %d of %d based on keep list",
            len(busco_runs),
            original_count,
        )

        if not busco_runs:
            LOGGER.error("No BUSCO runs matched keep list")
            return 1
# --------------------------------

    samples = sorted(busco_runs)
    LOGGER.info("Proceeding with %d BUSCO-labelled genomes", len(samples))

    if args.dry_run:
        LOGGER.info("DRY-RUN complete (BUSCO discovery only). Exiting before generating new files.")
        return 0

    gene_index = build_gene_index(busco_runs)
    if not gene_index:
        LOGGER.error("No single-copy BUSCO genes detected in provided runs")
        return 1

    total_samples = len(samples)
    fraction_threshold = ceil(total_samples * args.min_sample_fraction)
    min_required = max(args.min_samples, fraction_threshold)
    min_required = min(min_required, total_samples)
    LOGGER.info(
        "Per-gene minimum set to %d samples (%d total, %.0f%% coverage target)",
        min_required,
        total_samples,
        args.min_sample_fraction * 100,
    )

    per_gene_dir = args.out_dir / "01_single_copy_fastas"
    align_dir = args.out_dir / "02_alignments"
    concat_dir = args.out_dir / "03_concatenated"
    metadata_dir = args.out_dir / "metadata"
    tree_dir = args.out_dir / "04_trees"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    sample_stats_path = summarize_samples(gene_index, samples, metadata_dir, args.prefix)
    LOGGER.info("Sample stats saved to %s", sample_stats_path)

    # ============================
# SKIP 01 & 02, USE EXISTING ALIGNMENTS
# ============================

    if not align_dir.exists():
        LOGGER.error("02_alignments directory not found. Cannot proceed.")
        return 1

    LOGGER.info("Using existing alignments from %s", align_dir)

    aligned_fastas = sorted(align_dir.glob("*.aln.faa"))

    if not aligned_fastas:
        LOGGER.error("No alignment files found in 02_alignments")
        return 1


    concat_path, partition_path, gene_stats_path = concatenate_alignments(aligned_fastas, samples, concat_dir, args.prefix)
    LOGGER.info("Concatenated alignment: %s", concat_path)
    LOGGER.info("Partition file: %s", partition_path)
    LOGGER.info("Gene stats: %s", gene_stats_path)

    active_tree_programs = [prog for prog in args.tree_programs if prog != "none"]
    if not active_tree_programs:
        LOGGER.info("Skipping tree inference (no programs requested)")
        return 0

    for program in active_tree_programs:
        if program == "iqtree3":
            run_iqtree(concat_path, partition_path, tree_dir, args)
        elif program == "raxml-ng":
            run_raxml(concat_path, partition_path, tree_dir, args)
        else:
            LOGGER.warning("Unknown tree program %s", program)

    LOGGER.info("Pipeline completed successfully")
    return 0


if __name__ == "__main__":
    sys.exit(main())
EOF

