import argparse
import json
import re
import os
import pandas as pd

def read_and_process(directory):
    data = []
    for filename in os.listdir(directory):
        if filename.endswith('_timing.txt'):
            parameters = parse_filename(filename)
            with open(directory + '/' + filename, 'r') as f:
                timing = f.read()
                match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.*)', timing) 
                elapsed_time = match.groups()[0]
                elapsed_time = sum(x * float(t) for x, t in zip([1, 60, 3600], elapsed_time.split(":")[::-1]))

            res = parameters.copy()
            res['elapsed_time'] = elapsed_time
            if res['type'] == 'fastbe':
                fastbe_file = filename.replace('timing.txt', 'results.json')
                with open(directory + '/' + fastbe_file, 'r') as f:
                    content = json.load(f)
                    res['objective'] = content['objective']
            else:
                projection_file = filename.replace('timing.txt', 'results.txt')
                with open(directory + '/' + projection_file, 'r') as f:
                    obj = float(next(f))
                    res['objective'] = obj

            data.append(res)
    df = pd.DataFrame(data)
    return df
 
def parse_filename(filename):
    base_name = filename.split(f'_results.json')[0]
    params = base_name.split('_')
    parameters = {
        'm': int(params[0][1:]), 
        'n': int(params[1][1:]), 
        's': int(params[2][1:]), 
        'c': int(params[3][1:]), 
        'r': int(params[4][1:]),
        'type': params[5]
    }
    return parameters

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some json files into a pandas DataFrame.')
    parser.add_argument('directory', type=str, help='The directory where the json files are located.')
    parser.add_argument('--output', type=str, help='The output file name.')
    
    args = parser.parse_args()

    df = read_and_process(args.directory)
    if args.output:
        df.to_csv(args.output, index=False)


