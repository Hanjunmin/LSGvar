#!/usr/local/bin/python
# coding: iso-8859-1

import pandas as pd
import numpy as np
import multiprocessing
import sys
import re

def process_group(group_key, df, duplicates, output_prefix):
    group_df = df[df['group_key'] == group_key].copy() 

    try:
        extracted_data = group_df['col8'].str.extract(r'([^:]+):(\d+)-(\d+)')
        group_df[['chromosome', 'start', 'end']] = extracted_data
        group_df = group_df.dropna(subset=['start', 'end'])

        group_df['start'] = pd.to_numeric(group_df['start'], errors='coerce')
        group_df['end'] = pd.to_numeric(group_df['end'], errors='coerce')

        parts = group_key.split()
        if len(parts) == 3:
            chromosome, start, end = parts
            start = float(start)
            end = float(end)

            # dup
            if not duplicates[
                (duplicates['chromosome'] == chromosome) &
                (duplicates['start'] == start) &
                (duplicates['end'] == end)
            ].empty:
                # 过滤
                filtered_group_df = group_df[
                    (group_df['chromosome'].isin(duplicates['chromosome'])) &
                    (group_df['start'].isin(duplicates['start'])) &
                    (group_df['end'].isin(duplicates['end']))
                ]
                group_df = filtered_group_df
        else:
            return None  

        group_df = group_df.drop(columns=['chromosome', 'start', 'end'], errors='ignore')

    except Exception as e:
        print(f"Error processing group {group_key}: {e}")
        return None

    unique_samples = group_df['col8'].nunique()
    if unique_samples <= 1:
        return None

    frequencies = []
    for seq_name in group_df['col8'].unique():
        seq_df = group_df[group_df['col8'] == seq_name]
        min_pos = seq_df['col6'].min()
        max_pos = seq_df['col6'].max()
        length = max_pos - min_pos
        group_size = len(seq_df)
        frequency = group_size / length if length > 0 else 0.0
        frequencies.append((seq_name, frequency))

    frequencies_file = f"{output_prefix}_{group_key.replace(' ', '_')}_frequencies.txt"
    with open(frequencies_file, 'w') as f_freq:
        for seq_name, frequency in sorted(frequencies, key=lambda x: x[1], reverse=True):
            f_freq.write(f"{seq_name},{frequency:.10f}\n")

    frequencies_df = pd.read_csv(frequencies_file, sep=',', header=None, names=['sample', 'frequency'])
    del_file = f"{output_prefix}_{group_key.replace(' ', '_')}_del.txt"
    with open(del_file, 'w') as f_del:
        if (frequencies_df['frequency'] == 0).all():
            return None
        else:
            try:
                sample_del = frequencies_df.sort_values(by='frequency', ascending=True)['sample'].tolist()[1:]
                for sample in sample_del:
                    del_values = group_df[group_df['col8'] == sample]['col7'].tolist()
                    for val in del_values:
                        f_del.write(f"{val}\n")
            except Exception as e:
                print(f"Error processing group {group_key}: {e}")
                return None

    return group_key

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

    df = pd.read_csv(input_file, sep='\t', header=None,
                     names=['col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8'])

    df['group_key'] = df[['col1', 'col2', 'col3']].astype(str).agg(' '.join, axis=1)
    unique_group_keys = df['group_key'].unique()

    duplicates = pd.read_csv('duplicates.txt', sep='\t', header=None,
                             names=['chromosome', 'start', 'end'])
    duplicates['start'] = pd.to_numeric(duplicates['start'], errors='coerce')
    duplicates['end'] = pd.to_numeric(duplicates['end'], errors='coerce')

    #multiprocessing.Pool
    with multiprocessing.Pool(processes=60) as pool:
        results = pool.starmap(
            process_group,
            [(key, df, duplicates, output_prefix) for key in unique_group_keys]
        )
