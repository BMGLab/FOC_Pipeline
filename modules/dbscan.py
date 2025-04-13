import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import argparse

def read_gff(gff_file):
    """Reads a GFF file and extracts genomic coordinates."""
    data = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attributes = fields
            if feature == "gene":
                data.append([chrom, int(start), int(end)])
    
    return pd.DataFrame(data, columns=['chromosome', 'start', 'end'])

def perform_dbscan(df, eps=1000, min_samples=3):
    """Performs DBSCAN clustering on genomic data."""
    df['midpoint'] = (df['start'] + df['end']) / 2
    coordinates = np.array(df[['midpoint']])
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coordinates)
    df['cluster'] = clustering.labels_
    return df

def main():
    parser = argparse.ArgumentParser(description='DBSCAN clustering on gene annotations.')
    parser.add_argument('--input', required=True, help='Input GFF file')
    parser.add_argument('--output', required=True, help='Output clustered file')
    parser.add_argument('--eps', type=float, default=1000, help='DBSCAN epsilon value')
    parser.add_argument('--min_samples', type=int, default=3, help='Minimum samples for DBSCAN')
    args = parser.parse_args()

    df = read_gff(args.input)
    clustered_df = perform_dbscan(df, eps=args.eps, min_samples=args.min_samples)
    clustered_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
