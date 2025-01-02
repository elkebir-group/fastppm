import numpy as np
import pandas as pd
import networkx as nx
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Make input for Sapling')
    parser.add_argument('variant_matrix', type=str, help='Variant read matrix')
    parser.add_argument('total_matrix', type=str, help='Total read matrix')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    variant_matrix = np.loadtxt(args.variant_matrix).astype(int)
    total_matrix   = np.loadtxt(args.total_matrix).astype(int)

    print("sample_index\tmutation_index\tcluster_index\tvar\tdepth")
    for i in range(variant_matrix.shape[0]):
        for j in range(variant_matrix.shape[1]):
            print(f"{i}\t{j}\t{j}\t{variant_matrix[i,j]}\t{total_matrix[i,j]}")
    
