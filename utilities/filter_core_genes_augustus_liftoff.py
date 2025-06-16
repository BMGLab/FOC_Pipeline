import os
import re
import argparse
from math import ceil
from pathlib import Path
from collections import defaultdict


def extract_gene_names_from_faa(faa_path):
    genes = set()
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                genes.add(line.split()[0][1:])
    return genes

def extract_gene_ids_from_gff(gff_path):
    genes = set()
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or '\tgene\t' not in line:
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            match = re.search(r'ID=([^;]+)', fields[8])
            if match:
                genes.add(match.group(1))
    return genes

def count_core_genes(sample_dirs, use_faa):
    gene_occurrence = defaultdict(int)
    sample_count = 0

    for path in sample_dirs:
        gff = next(path.glob("*.gff"), None)
        if not gff:
            continue

        gff_genes = extract_gene_ids_from_gff(gff)

        if use_faa:
            faa = next(path.glob("*.faa"), None)
            if not faa:
                print(f"⚠️ Missing .faa for {path.name}")
                continue
            faa_genes = extract_gene_names_from_faa(faa)
            shared = gff_genes & faa_genes
        else:
            shared = gff_genes

        for gene in shared:
            gene_occurrence[gene] += 1

        sample_count += 1

    return gene_occurrence, sample_count


def write_core_gene_list(gene_occurrence, sample_count, threshold=0.7):
    min_required = ceil(threshold * sample_count)
    return set(g for g, c in gene_occurrence.items() if c >= min_required)


def filter_faa(input_path, output_path, gene_set):
    with open(input_path) as fin, open(output_path, 'w') as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                gene = line.split()[0][1:]
                keep = gene in gene_set
            if keep:
                fout.write(line)

def filter_gff(input_path, output_path, gene_set):
    with open(input_path) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith("#") or '\tgene\t' not in line:
                continue
            match = re.search(r'ID=([^;]+)', line)
            if match and match.group(1).rstrip() in gene_set:
                fout.write(line)

def run_filtering(root_dir, use_faa=True, output_dir="core_filtered"):
    root = Path(root_dir)
    samples = [p for p in root.iterdir() if p.is_dir()]
    print(f"🔍 Found {len(samples)} samples in {root_dir}")

    gene_occurrence, sample_count = count_core_genes(samples, use_faa)
    if sample_count == 0:
        print("❌ No valid samples found.")
        return

    core_genes = write_core_gene_list(gene_occurrence, sample_count)
    os.makedirs(output_dir, exist_ok=True)

    for path in samples:
        name = path.name
        gff_in = next(path.glob("*.gff"), None)
        faa_in = next(path.glob("*.faa"), None) if use_faa else None

        if gff_in:
            gff_out = Path(output_dir) / f"{name}_core.gff"
            filter_gff(gff_in, gff_out, core_genes)

        if use_faa and faa_in:
            faa_out = Path(output_dir) / f"{name}_core.faa"
            filter_faa(faa_in, faa_out, core_genes)

    print(f"✅ Filtered outputs written to: {output_dir}")
    print(f"🧬 Core gene count: {len(core_genes)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter core genes from Liftoff or AUGUSTUS output.")
    parser.add_argument("input_path", help="Path to either data/foc/augustus or data/foc/liftoff")
    parser.add_argument("--output_dir", default="core_filtered", help="Output directory")
    parser.add_argument("--no_faa", action="store_true", help="Set this if only GFF files exist (e.g., Liftoff)")

    args = parser.parse_args()
    run_filtering(args.input_path, use_faa=not args.no_faa, output_dir=args.output_dir)