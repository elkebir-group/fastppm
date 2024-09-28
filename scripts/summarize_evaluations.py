import re
import os
import argparse
import networkx as nx
import pandas as pd

def load_files(directory):
    data = []
    for algorithm in os.listdir(directory):
        for subdir in os.listdir(os.path.join(directory, algorithm)):
            match = re.search(r'n(\d+)_s(\d+)_c(\d+)_r(\d+)', subdir)
            n, s, c, r = match.groups()
            with open(os.path.join(directory, algorithm, subdir, 'timing.txt')) as f:
                timing = f.read().strip()
                match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', timing)
                elapsed_time = match.groups()[0]
                elapsed_time = sum(x * float(t) for x, t in zip([1, 60, 3600], elapsed_time.split(":")[::-1]))

            true_tree     = nx.read_adjlist(os.path.join('data/simulations/', subdir, 'sim_tree.txt'))
            inferred_tree = nx.read_adjlist(os.path.join(directory, algorithm, subdir, 'tree.txt'))

            true_positives = len(set(true_tree.edges()) & set(inferred_tree.edges()))
            false_positives = len(set(inferred_tree.edges()) - set(true_tree.edges()))
            false_negatives = len(set(true_tree.edges()) - set(inferred_tree.edges()))
            f1_score = 2 * true_positives / (2 * true_positives + false_positives + false_negatives)

            data.append({
                'algorithm': algorithm,
                'n': int(n),
                's': int(s),
                'c': int(c),
                'r': int(r),
                'elapsed_time': elapsed_time,
                'true_positives': true_positives,
                'false_positives': false_positives,
                'false_negatives': false_negatives,
                'f1_score': f1_score,
            })
            
    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Process the results of the evaluations.')
    parser.add_argument('directory', type=str, help='The directory containing the results.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('summary.csv', index=False)

if __name__ == '__main__':
    main()

