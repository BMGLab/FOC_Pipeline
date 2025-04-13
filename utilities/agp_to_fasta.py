#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_agp(agp_file):
    """Parse AGP file and return scaffold structure"""
    scaffolds = {}
    with open(agp_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if parts[4] == 'W':  # Only process sequence lines
                scaffold_name = parts[0]
                component_name = parts[5]
                strand = parts[8]
                if scaffold_name not in scaffolds:
                    scaffolds[scaffold_name] = []
                scaffolds[scaffold_name].append((component_name, strand))
    return scaffolds

def convert_agp_to_fasta(agp_file, input_fasta, output_fasta):
    """Convert AGP file to FASTA format"""
    # Read input sequences
    sequences = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    
    # Parse AGP file
    scaffolds = parse_agp(agp_file)
    
    # Generate new sequences
    new_records = []
    for scaffold_name, components in scaffolds.items():
        scaffold_seq = ""
        for comp_name, strand in components:
            if comp_name in sequences:
                seq = str(sequences[comp_name].seq)
                if strand == '-':
                    seq = str(Seq(seq).reverse_complement())
                scaffold_seq += seq
        
        record = SeqRecord(
            Seq(scaffold_seq),
            id=scaffold_name,
            description=""
        )
        new_records.append(record)
    
    # Write output FASTA
    SeqIO.write(new_records, output_fasta, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert AGP file to FASTA format")
    parser.add_argument("agp_file", help="Input AGP file")
    parser.add_argument("input_fasta", help="Input FASTA file with original sequences")
    parser.add_argument("output_fasta", help="Output FASTA file")
    args = parser.parse_args()
    
    convert_agp_to_fasta(args.agp_file, args.input_fasta, args.output_fasta)