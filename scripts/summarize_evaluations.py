import re
import os
import argparse
import numpy as np
import networkx as nx
import pandas as pd
import fastppm

SIMULATION_DIR = "/n/fs/ragr-research/projects/fastppm-data/simulations"

def load_files(directory):
    data = []
    for algorithm in os.listdir(directory):
        for subdir in os.listdir(os.path.join(directory, algorithm)):
            match = re.search(r'n(\d+)_s(\d+)_c(\d+)_r(\d+)', subdir)
            n, s, c, r = match.groups()

            if not os.path.exists(os.path.join(directory, algorithm, subdir, 'timing.txt')):
                continue

            with open(os.path.join(directory, algorithm, subdir, 'timing.txt')) as f:
                timing = f.read().strip()
                match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', timing)
                elapsed_time = match.groups()[0]
                elapsed_time = sum(x * float(t) for x, t in zip([1, 60, 3600], elapsed_time.split(":")[::-1]))

            true_tree = nx.read_adjlist(os.path.join(SIMULATION_DIR, subdir, 'sim_tree.txt'), create_using=nx.DiGraph())

            if 'sapling' in algorithm:
                inferred_tree_fname = os.path.join(directory, algorithm, subdir, 'sapling_output.txt')
                edges = []
                with open(inferred_tree_fname) as f:
                    next(f)
                    next(f)
                    for line in f:
                        if line.strip() == '':
                            break
                        u, v = line.strip().split()
                        if u == '-1' or v == '-1':
                            continue
                        edges.append((u, v))
                inferred_tree = nx.DiGraph(edges)
            else:
                inferred_tree = nx.read_adjlist(os.path.join(directory, algorithm, subdir, 'tree.txt'), create_using=nx.DiGraph())

            # results_file  = os.path.join(directory, algorithm, subdir, 'results.npz')
            # results_npz   = np.load(results_file)
            # log_likelihood_orchard = results_npz['llh'][0]

            valid_inferred_tree = True
            if inferred_tree.number_of_nodes() != true_tree.number_of_nodes():
                valid_inferred_tree = False

            if inferred_tree.number_of_edges() != true_tree.number_of_edges():
                valid_inferred_tree = False

            negative_log_likelihood = None
            if valid_inferred_tree:
                variant_matrix = np.loadtxt(os.path.join(SIMULATION_DIR, subdir, 'sim_variant_matrix.txt')).astype(int)
                total_matrix   = np.loadtxt(os.path.join(SIMULATION_DIR, subdir, 'sim_total_matrix.txt')).astype(int)

                n = variant_matrix.shape[1]
                adj_list = [[] for _ in range(n)]
                for (i, j) in inferred_tree.edges():
                    adj_list[int(i)].append(int(j))
    
                try:
                    res = fastppm.regress_counts(adj_list, variant_matrix.tolist(), total_matrix.tolist(), loss_function='binomial')
                    negative_log_likelihood = 0
                    #negative_log_likelihood = res['objective']
                except:
                    print('Error in', algorithm, subdir)
                    valid_inferred_tree = False

            true_positives = len(set(true_tree.edges()) & set(inferred_tree.edges()))
            false_positives = len(set(inferred_tree.edges()) - set(true_tree.edges()))
            false_negatives = len(set(true_tree.edges()) - set(inferred_tree.edges()))
            f1_score = 2 * true_positives / (2 * true_positives + false_positives + false_negatives)

            true_tree_closure = nx.transitive_closure(true_tree)
            inferred_tree_closure = nx.transitive_closure(inferred_tree)
            ancestor_descendant_distance = nx.symmetric_difference(true_tree_closure, inferred_tree_closure).number_of_edges()
            normalizing_factor = len(set(true_tree_closure.edges()) | set(inferred_tree_closure.edges()))
            normalized_ancestor_descendant_distance = ancestor_descendant_distance / normalizing_factor

            closure_true_positives = len(set(true_tree_closure.edges()) & set(inferred_tree_closure.edges()))
            closure_false_positives = len(set(inferred_tree_closure.edges()) - set(true_tree_closure.edges()))
            closure_false_negatives = len(set(true_tree_closure.edges()) - set(inferred_tree_closure.edges()))
            closure_f1_score = 2 * closure_true_positives / (2 * closure_true_positives + closure_false_positives + closure_false_negatives)

            row = {
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
                'valid_inferred_tree': valid_inferred_tree,
                'negative_log_likelihood': negative_log_likelihood,
                'ancestor_descendant_distance': ancestor_descendant_distance,
                'normalized_ancestor_descendant_distance': normalized_ancestor_descendant_distance,
                'f1_score_closure': closure_f1_score,
            }

            print(row)
            data.append(row)
            
    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description='Process the results of the evaluations.')
    parser.add_argument('directory', type=str, help='The directory containing the results.')
    args = parser.parse_args()

    df = load_files(args.directory)
    df.to_csv('summary.csv', index=False)

if __name__ == '__main__':
    main()

