#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  chromosome_similarity.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>


import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import parse_pangenome


def organize_data(parsed_data):
    chromosome_data = defaultdict(lambda: defaultdict(list))
    for node, data in parsed_data.items():
        chromosome = data['chromosome']
        node_length = int(data['node_length'])
        samples = data['samples']
        for sample in samples:
            chromosome_data[chromosome][sample].append((node, node_length))
    return chromosome_data


def calculate_windows(nodes, window_size):
    windowed_data = []
    current_window = []
    current_length = 0
    window_start = 0
    for node, length in nodes:
        while current_length + length > window_size:
            windowed_data.append((window_start, current_length, current_window))
            window_start += window_size
            current_length -= window_size
            current_window = [(n, l) for n, l in current_window if current_length + l > window_start]
        current_window.append((node, length))
        current_length += length
    if current_length > 0:
        windowed_data.append((window_start, current_length, current_window))

    return windowed_data


def determine_similarity(sample, windows):
    similarity_data = []
    for start, length, window_nodes in windows:
        sample_counts = defaultdict(int)
        for node, node_length in window_nodes:
            node_samples = parsed_data[node]['samples']
            for query in node_samples:
                if query.split('.')[0] != sample.split('.')[0] and query in haplotype_colors:
                    sample_counts[query] += node_length
        if sample_counts:
            most_similar_sample = max(sample_counts, key=sample_counts.get)
        else:
            most_similar_sample = 'unique'
        similarity_data.append((length, most_similar_sample))

    return pd.DataFrame(similarity_data, columns=['size', 'haplotype'])


def plot_haplotypes(output_name):
    fig, ax = plt.subplots(figsize=(10, 7))
    # Vertical space between each sample
    y_offset = 1
    size_list = []

    for sample in reversed(sample_order):
        df = output_dict[sample]
        size_list.append(df['size'].sum())
        start = 0  # Starting position for the first node in each sample
        for _, row in df.iterrows():
            haplotype = row['haplotype']
            size = row['size']
            color = haplotype_colors.get(haplotype, 'grey')
            # Draw a rectangle for each node
            ax.add_patch(plt.Rectangle((start, y_offset), size, 0.8, color=color))
            start += size
        # Label for the sample
        ax.text(-5000, y_offset + 0.4, sample, va='center', ha='right', fontsize=9)
        # Move to the next sample position
        y_offset += 1

    # Set limits and labels
    ax.set_xlim(0, max(size_list))
    ax.set_ylim(0.2, y_offset)
    ax.set_yticks([])
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xlabel('Chromosome Size (Mb)', fontsize=10)
    ax.set_title(output_name, fontsize=12)

    # Add a legend with a fixed location outside the figure
    haplotypes_in_legend = set()
    for haplotype, color in haplotype_colors.items():
        haplotypes_in_legend.add(haplotype)
        ax.plot([], [], color=color, label=haplotype)
    ax.legend(title='Best Match', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9, title_fontsize=10)

    plt.savefig(f'{output_dir}/{output_name}.pdf', format='pdf', bbox_inches='tight', dpi=400)
    plt.close()


def parse_list_file(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]


def parse_color_tsv(file_path):
    colors = {}
    with open(file_path) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    raise argparse.ArgumentTypeError(f"Invalid line in {file_path}: {line}")
                hap, color = parts
                colors[hap] = color
    colors.setdefault('unique', 'grey')
    return colors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate chromosome similarity heatmap.')

    parser.add_argument('--input_table', type=str, required=True,
                        help='Path to the input TSV file with node information.')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Directory to save output files.')
    parser.add_argument('--chrom', type=str, required=True,
                        help='Chromosome to analyze.')
    parser.add_argument('--window_size', type=int, default=200000,
                        help='Window size for aggregation (default: 200000).')

    parser.add_argument('--haplotype_colors_file', type=parse_color_tsv, required=True,
                        help='TSV file with haplotype and hex color code.')
    parser.add_argument('--sample_order_file', type=parse_list_file, required=True,
                        help='File with sample drawing order, one sample per line.')

    args = parser.parse_args()

    input_table = args.input_table
    output_dir = args.output_dir
    chrom = args.chrom
    window_size = args.window_size
    haplotype_colors = args.haplotype_colors_file
    sample_order = args.sample_order_file

    if not os.path.exists(os.path.dirname(output_dir + '/')):
        os.makedirs(os.path.dirname(output_dir + '/'))

    parsed_data, sample_list, _ = parse_pangenome.parse_table(input_table, chrom)
    chromosome_data = organize_data(parsed_data)

    for chromosome, samples in chromosome_data.items():
        output_dict = {}
        for sample, nodes in samples.items():
            print(sample)
            windowed_data = calculate_windows(nodes, window_size)
            similarity_df = determine_similarity(sample, windowed_data)
            output_dict[sample] = similarity_df

        plot_haplotypes(f'{chromosome}_similarity_heatmap_{window_size}')
