import re
import json
import os
import argparse
import numpy as np
import networkx as nx
import pandas as pd

simulated_data_dir = "/n/fs/ragr-research/projects/fastppm-data/simulations"

def process_instance(directory, algorithm, subdir):
    if len(subdir.split('_')) == 4:
        n, s, c, r = subdir.split('_')
        n, s, c, r = int(n[1:]), int(s[1:]), int(c[1:]), int(r[1:])
        k = None

    elif len(subdir.split('_')) == 5:
        n, s, c, r, k = subdir.split('_')
        n, s, c, r, k = int(n[1:]), int(s[1:]), int(c[1:]), int(r[1:]), int(k[1:])

    if not os.path.exists(os.path.join(directory, algorithm, subdir, 'timing.txt')):
        return None 

    with open(os.path.join(directory, algorithm, subdir, 'timing.txt')) as f:
        timing = f.read().strip()
        match = re.search(r'User time \(seconds\): ([0-9]+\.[0-9]+)', timing)
        user_time = float(match.group(1))
        match = re.search(r'System time \(seconds\): ([0-9]+\.[0-9]+)', timing)
        system_time = float(match.group(1))
        elapsed_time = user_time + system_time

    nll = None
    if algorithm == "projection_l2":
        with open(os.path.join(directory, algorithm, subdir, 'output.txt')) as f:
            objective = float(next(f).strip())
    else:
        with open(os.path.join(directory, algorithm, subdir, 'output.json')) as f:
            print(os.path.join(directory, algorithm, subdir, 'output.json'))
            res = json.load(f)
            objective = res['objective']

        if 'l2' not in algorithm:
            frequency_matrix = np.array(res['frequency_matrix'])
            variant_matrix   = np.loadtxt(f"{simulated_data_dir}/n{n}_s{s}_c{c}_r{r}/sim_variant_matrix.txt")
            total_matrix     = np.loadtxt(f"{simulated_data_dir}/n{n}_s{s}_c{c}_r{r}/sim_total_matrix.txt")
            nll = np.log(frequency_matrix, out=np.zeros_like(frequency_matrix), where=(frequency_matrix!=0)) * variant_matrix 
            nll += (total_matrix - variant_matrix) * np.log(1 - frequency_matrix, out=np.zeros_like(frequency_matrix), where=(frequency_matrix!=1))
            nll = -np.sum(nll)
            print(nll)

    return {
        'n': n,
        's': s,
        'c': c,
        'r': r,
        'k': k,
        'algorithm': algorithm,
        'elapsed_time': elapsed_time,
        'objective': objective,
        'nll': nll
    }

def main():
    parser = argparse.ArgumentParser(description='Process the results of the evaluations.')
    parser.add_argument('directory', type=str, help='The directory containing the results.')
    args = parser.parse_args()

    rows = []
    for algorithm in os.listdir(args.directory):
        for subdir in os.listdir(os.path.join(args.directory, algorithm)):
            row = process_instance(args.directory, algorithm, subdir)
            if row is None: continue
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv('summary_regressions.csv', index=False)

if __name__ == '__main__':
    main()

