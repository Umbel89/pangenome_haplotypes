#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  chromosome_structure.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import os
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import parse_pangenome


def append_dictionary():
    global chrom_dict, node_dict, matrix_dict
    chrom_dict = {}
    node_dict = {}
    matrix_dict = {}

    for chromosome in chrom_list:
        chrom_dict[chromosome] = {}
        node_dict[chromosome] = {}
        matrix_dict[chromosome] = {}
        for sample in sample_list:
            chrom_dict[chromosome][sample] = {}
            matrix_dict[chromosome][sample] = []


def graph_variation():
    for node in parsed_data:
        node_length = int(parsed_data[node]['node_length'])
        chromosome = parsed_data[node]['chromosome']
        node_samples = parsed_data[node]['samples']
        haplotype = assign_haplotype(node_samples)
        node_dict[chromosome][node] = haplotype
        for sample in sample_list:
            if sample not in node_samples:
                sample_haplotype = 'absent'
                matrix_value = 0
            else:
                sample_haplotype = haplotype
                matrix_value = node_length
            matrix_dict[chromosome][sample].append(matrix_value)
            chrom_dict[chromosome][sample][node] = {'size': node_length, 'haplotype': sample_haplotype}


def assign_haplotype(node_samples):
    if len(node_samples) == len(sample_list):
        return "core"
    elif len(node_samples) == 1:
        return "unique"
    # check if both samples in the list are the two phases of one sample
    elif len(node_samples) == 2:
        if node_samples[0].split('.')[0] == node_samples[1].split('.')[0]:
            return "unique"

    for sample in haplotype_order:
        if sample in node_samples:
            return sample


def calculate_windows(input_dict):
    variation_df = pd.DataFrame(input_dict).T
    variation_df.sort_index(inplace=True)
    results = []
    current_window_start = 0
    current_window_end = window_size
    current_window_size = 0
    haplotype_counts = {}

    for _, row in variation_df.iterrows():
        size, haplotype = row['size'], row['haplotype']
        while size > 0:
            remaining_space = window_size - current_window_size
            if size <= remaining_space:
                add_size = size
            else:
                add_size = remaining_space
            current_window_size += add_size
            size -= add_size
            if haplotype not in haplotype_counts:
                haplotype_counts[haplotype] = 0
            haplotype_counts[haplotype] += add_size

            if current_window_size == window_size:
                find_best_match(current_window_size, haplotype_counts, results)
                current_window_start += window_size
                current_window_end += window_size
                current_window_size = 0
                haplotype_counts = {}

    # if there is a smaller window left at the end of the chromosome
    if current_window_size > 0:
        find_best_match(current_window_size, haplotype_counts, results)

    return pd.DataFrame(results, columns=['size', 'haplotype'])


def find_best_match(current_window_size, haplotype_counts, results):
    # filter out core, unless only core nodes are present in this window
    if len(haplotype_counts) == 1 and "core" in haplotype_counts:
        max_haplotype = "core"
    else:
        filtered_haplotypes = {k: v for k, v in haplotype_counts.items() if k != "core"}
        max_haplotype, max_count = max(filtered_haplotypes.items(), key=lambda item: item[1])
        if max_count < args.threshold*window_size:
            max_haplotype = "core"

    results.append({'size': current_window_size, 'haplotype': max_haplotype})


def plot_haplotypes(output_name):
    fig, ax = plt.subplots(figsize=(10, 6))
    # Vertical space between each sample
    y_offset = 1

    for sample in reversed(heatmap_order):
        df = output_dict[sample]
        start = 0  # Starting position for the first node in each sample
        for _, row in df.iterrows():
            haplotype = row['haplotype']
            size = row['size'] #/ total_size
            color = haplotype_colors.get(haplotype, 'grey')
            # Draw a rectangle for each node
            ax.add_patch(plt.Rectangle((start, y_offset), size, 0.8, color=color))
            start += size
        # Label for the sample
        ax.text(-5000, y_offset + 0.4, sample, va='center', ha='right', fontsize=9)
        # Move to the next sample position
        y_offset += 1

    # Set limits and labels
    ax.set_xlim(0, df['size'].sum()+0.01*df['size'].sum())
    ax.set_ylim(0.2, y_offset)
    ax.set_yticks([])
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xlabel('Pangenome Size (MB)', fontsize=10)
    ax.set_title(output_name, fontsize=12)

    # Add a legend with a fixed location outside the figure
    haplotypes_in_legend = set()
    for haplotype, color in haplotype_colors.items():
        haplotypes_in_legend.add(haplotype)
        ax.plot([], [], color=color, label=haplotype)
    ax.legend(title='Haplotype', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9, title_fontsize=10)

    plt.savefig(f'{args.output_dir}/{output_name}.pdf', format='pdf', bbox_inches='tight', dpi=400)
    plt.close()


def write_dataframe(input_df, output_name):
    # write dataframe to file
    output_fn = f'{args.output_dir}/{output_name}.tsv'
    input_df.to_csv(output_fn, sep='\t')


def write_csv(output_name):
    output_fn = f'{args.output_dir}/{output_name}.csv'
    if not os.path.isfile(output_fn) or os.path.getsize(output_fn) == 0:
        with open(output_fn, 'w') as the_file:
            the_file.write('Name,Chromosome,Race,Colour\n')
            for chromosome in node_dict:
                first_node = True
                for node, haplotype in node_dict[chromosome].items():
                    colour = haplotype_colors.get(haplotype, 'grey')
                    node = str(int(node) - 1139705)
                    if first_node:
                        the_file.write(f'{node},{chromosome},{haplotype},{colour}\n')
                        first_node = False
                    else:
                        the_file.write(f'{node},,{haplotype},{colour}\n')


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
    # Ensure defaults
    colors.setdefault('core', '#43A2CA')
    colors.setdefault('unique', '#FFC000')
    colors.setdefault('absent', 'white')
    colors.setdefault('other', 'grey')
    return colors


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze chromosome structure and generate haplotype plots.')

    parser.add_argument('--input_table', type=str, required=True, help='Input table file path.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save output.')
    parser.add_argument('--chrom', type=str, required=True, help='Chromosome to analyze.')

    parser.add_argument('--window_size', type=int, default=4000, help='Size of analysis window.')
    parser.add_argument('--threshold', type=float, default=0.05, help='Threshold for haplotype assignment.')
    parser.add_argument('--min_node', type=int, default=0, help='Minimum node to include.')
    parser.add_argument('--max_node', type=float, default=float('inf'), help='Maximum node to include.')

    parser.add_argument('--haplotype_order_file', type=parse_list_file, required=True,
                        help='File with haplotype order, one per line.')
    parser.add_argument('--haplotype_colors_file', type=parse_color_tsv, required=True,
                        help='TSV file with haplotype and color code.')
    parser.add_argument('--heatmap_order_file', type=parse_list_file, required=True,
                        help='File with heatmap sample order, one per line.')


    args = parser.parse_args()

    window_size = args.window_size
    if not os.path.exists(os.path.dirname(args.output_dir + '/')):
        os.makedirs(os.path.dirname(args.output_dir + '/'))

    parsed_data, sample_list, chrom_list = parse_pangenome.parse_table(args.input_table, args.chrom, 
                                                                       args.min_node, args.max_node)

    output_dict = {}
    append_dictionary()
    graph_variation()
    write_csv(f'{args.chrom}_node_haplotypes_{args.threshold}')

    for chromosome, sample_dict in chrom_dict.items():
        print(chromosome)
        output_dict = {}
        for sample, location_dict in sample_dict.items():
            print(sample)
            windows_df = calculate_windows(location_dict)
            output_dict[sample] = windows_df
        plot_haplotypes(f'{chromosome}_accessory_haplotypes_{window_size}_{args.threshold}')
