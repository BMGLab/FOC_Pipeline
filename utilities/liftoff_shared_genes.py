# 🔧 What it does:
# 	1.	Parses multiple .gff files to extract gene names (Name= attribute).
# 	2.	Counts how many times each gene name appears across files.
# 	3.	Filters genes that are present in at least 8 files.
# 	4.	Extracts full GFF records for those shared genes, writing per-sample filtered GFF files.

# ✅ How to run it:
# 	1.	Save the code above as liftoff_shared_genes.py.
# 	2.	Ensure you have your *.gff files under ~/liftoff/*/*predicted_genes.gff.
# 	3.	Run from terminal: python3 liftoff_shared_genes.py

import os
import glob
from collections import defaultdict

# Path to Liftoff output files
input_dir = os.path.expanduser('~/liftoff/')
gff_files = glob.glob(os.path.join(input_dir, '*/*predicted_genes.gff'))

# Step 1: Count how many files each gene name appears in
gene_file_count = defaultdict(set)  # gene_name -> set of filenames

for filepath in gff_files:
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or '\tgene\t' not in line:
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            for attr in attributes.split(';'):
                if attr.strip().startswith('Name='):
                    gene_name = attr.split('=')[1].strip()
                    gene_file_count[gene_name].add(filepath)
                    break  # One Name= per gene assumed

# Filter genes shared in 8 or more isolates
shared_genes = {gene for gene, files in gene_file_count.items() if len(files) >= 8}
print(f"Found {len(shared_genes)} shared genes in 8 or more isolates.")

# Step 2: Filter each GFF file to keep only shared genes
for filepath in gff_files:
    sample_name = os.path.basename(filepath).replace('_predicted_genes.gff', '')
    output_path = os.path.join(os.path.dirname(filepath), f"{sample_name}_shared_8plus.gff")

    with open(filepath) as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            if '\tgene\t' not in line:
                continue
            fields = line.strip().split('\t')
            attributes = fields[8]
            keep = False
            for attr in attributes.split(';'):
                if attr.strip().startswith('Name='):
                    gene_name = attr.split('=')[1].strip()
                    if gene_name in shared_genes:
                        keep = True
                    break
            if keep:
                fout.write(line)

print("Filtered GFF files written.")