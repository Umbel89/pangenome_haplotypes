#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_pangenome.py
#
#  Copyright 2021 Petros Skiadas <p.skiadas@uu.nl>

import csv
import ast
from collections import OrderedDict


def parse_table(file_path, chrom, min_node=0, max_node=float("inf")):
    data = {}
    chrom_list = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)  # Read the header
        headers[0] = 'node_id'  # Give a name to the first column header
        sample_list = headers[8:]

        for row in reader:
            node_number = int(row[0])  # The first column is used as the key
            if min_node <= node_number <= max_node:
                row_data = dict(zip(headers[1:], row[1:]))  # Create a dictionary for the rest of the columns
                node_chrom = row_data['chromosome']
                # check if the chromosome matches the input chrom
                if node_chrom == chrom or chrom == 'all':
                    if node_chrom not in chrom_list:
                        chrom_list.append(node_chrom)
                    # Convert specific fields from string representation to Python objects
                    row_data['repeat'] = ast.literal_eval(row_data['repeat'])
                    row_data['gene_id'] = ast.literal_eval(row_data['gene_id'])
                    row_data['samples'] = row_data['samples'].split(',')

                    # Convert all Pe fields
                    for field in row_data.keys():
                        if field.startswith('Pe') and row_data[field]:
                            row_data[field] = ast.literal_eval(row_data[field])
                        elif field.startswith('Pe') and not row_data[field]:
                            row_data[field] = []

                    data[node_number] = row_data

    data_ordered = OrderedDict(sorted(data.items()))

    return data_ordered, sample_list, chrom_list
