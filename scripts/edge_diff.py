import argparse

def read_edge_set(filename):
    edge_set = set()
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split()
            if not parts:
                continue  # Skip empty lines
            u = int(parts[0])
            for v_str in parts[1:]:
                v = int(v_str)
                edge_set.add((u, v))
    return edge_set

def main():
    parser = argparse.ArgumentParser(description='Compute edge set symmetric differences.')
    parser.add_argument('files', nargs='+', help='Input files in adjacency list format.')
    args = parser.parse_args()

    files = args.files
    if len(files) < 2:
        print("Please provide at least two files.")
        return

    # Read the reference edge set
    reference_file = files[0]
    ref_edge_set = read_edge_set(reference_file)
    print(f'Reference file: {reference_file}')
    print(f'Number of edges in reference: {len(ref_edge_set)}')

    # For each other file, compare to reference
    for filename in files[1:]:
        edge_set = read_edge_set(filename)
        print(f'\nComparing to file: {filename}')
        print(f'Number of edges: {len(edge_set)}')

        # Compute symmetric difference
        sym_diff = ref_edge_set.symmetric_difference(edge_set)
        num_sym_diff = len(sym_diff)
        print(f'Number of edges in symmetric difference: {num_sym_diff}')

        # False positives: edges in edge_set but not in ref_edge_set
        false_positives = edge_set - ref_edge_set
        num_false_positives = len(false_positives)
        print(f'Number of false positives: {num_false_positives}')

        # False negatives: edges in ref_edge_set but not in edge_set
        false_negatives = ref_edge_set - edge_set
        num_false_negatives = len(false_negatives)
        print(f'Number of false negatives: {num_false_negatives}')

if __name__ == '__main__':
    main()

