# 🔍 What the script does:
# 	1.	Reads multiple .faa (protein) files.
# 	2.	Extracts gene IDs (headers in FASTA format).
# 	3.	Identifies gene names that appear in 8 or more isolates.
# 	4.	Filters each .faa file to include only these shared genes.

# ✅ How to run it:
# 	1.	Save this code to a file: augustus_shared_genes.py
# 	2.	Make sure your .faa files are in folders like:
#       ~/augustus/10/10_augustus_predictions.faa, etc.
# 	3.	Then run:
#       python3 augustus_shared_genes.py

import os
import glob
from collections import defaultdict

# Define the isolate IDs used
isolate_ids = [10, 13, 2, 34, 36, 37, 39, 41, 42, 43, 45]
input_dir = os.path.expanduser('~/augustus/')

# Step 1: Count gene appearances across isolates
gene_to_isolates = defaultdict(set)

for isolate_id in isolate_ids:
    faa_path = os.path.join(input_dir, str(isolate_id), f"{isolate_id}_augustus_predictions.faa")
    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                gene_name = line.split()[0][1:]  # Remove ">" and take first word
                gene_to_isolates[gene_name].add(isolate_id)

# Filter genes that appear in 8 or more isolates
shared_genes = {gene for gene, ids in gene_to_isolates.items() if len(ids) >= 8}
print(f"Found {len(shared_genes)} shared genes in 8 or more isolates.")

# Step 2: Filter each .faa file for shared genes
for isolate_id in isolate_ids:
    input_faa = os.path.join(input_dir, str(isolate_id), f"{isolate_id}_augustus_predictions.faa")
    output_faa = os.path.join(input_dir, str(isolate_id), f"{isolate_id}_shared_8plus.faa")

    with open(input_faa) as fin, open(output_faa, 'w') as fout:
        keep = False
        for line in fin:
            if line.startswith('>'):
                gene_name = line.split()[0][1:]
                keep = gene_name in shared_genes
                if keep:
                    fout.write(line)
            else:
                if keep:
                    fout.write(line)

print("Filtered FAA files written.")