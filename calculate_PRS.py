#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd

def calculate_PRS(snps_weights_file, genotypes_file, out_file):
    # Load SNPs and weights
    snps_weights = pd.read_csv(snps_weights_file, sep='\t', header=None, names=['SNP', 'Weight'])
    
    # Load genotypes
    genotypes = pd.read_csv(genotypes_file + '.bed', sep='\t', header=None, dtype=np.uint8)
    
    # Convert genotypes to 0, 1, 2 coding
    genotypes = genotypes.iloc[:, 4:].values.reshape(-1, 2)
    genotypes = np.sum(genotypes, axis=1)
    
    # Calculate polygenic risk score
    PRS = np.dot(snps_weights['Weight'], genotypes)
    
    # Save polygenic risk score to file
    with open(out_file, 'w') as f:
        f.write(str(PRS))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate polygenic risk score')
    parser.add_argument('--snps_weights', type=str, help='SNPs and weights file')
    parser.add_argument('--genotypes', type=str, help='Genotypes file prefix')
    parser.add_argument('--out', type=str, help='Output file')
    args = parser.parse_args()

    calculate_PRS(args.snps_weights, args.genotypes, args.out)
