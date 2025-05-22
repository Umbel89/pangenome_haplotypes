# Pangenome Haplotypes
Code to visualise the haplotypes of a region of a pangenome graph. The scripts take as input a table of a pangenome graph from the pipeline https://github.com/TeamMGE/Skiadas2024_pangenome_archive, parses a given region of the table, calculates the haplotypes in the graph with a given priority, and plots the resulted haplotypes in a multiple sequence alignment.

![Pangenome Haplotypes](https://github.com/user-attachments/assets/ea759aaf-c7bb-4d38-9c6a-57b6f0bb6fb4)

## Installation
```bash
git clone https://github.com/Umbel89/pangenome_haplotypes.git
```

## Usage - Haplotype figure
The haplotype figure can be created by running the chromosome_structure.py script:
```
python "/mnt/WORKSPACE/psk_workspace/scripts/pangenome/chromosome_structure.py" -h
usage: chromosome_structure.py [-h] --input_table INPUT_TABLE --output_dir OUTPUT_DIR --chrom
                               CHROM [--window_size WINDOW_SIZE] [--threshold THRESHOLD]
                               [--min_node MIN_NODE] [--max_node MAX_NODE]
                               --haplotype_order_file HAPLOTYPE_ORDER_FILE
                               --haplotype_colors_file HAPLOTYPE_COLORS_FILE
                               --heatmap_order_file HEATMAP_ORDER_FILE

Analyze chromosome structure and generate haplotype plots.

options:
  -h, --help            show this help message and exit
  --input_table INPUT_TABLE
                        Input table file path.
  --output_dir OUTPUT_DIR
                        Directory to save output.
  --chrom CHROM         Chromosome to analyze.
  --window_size WINDOW_SIZE
                        Size of analysis window.
  --threshold THRESHOLD
                        Threshold for haplotype assignment.
  --min_node MIN_NODE   Minimum node to include.
  --max_node MAX_NODE   Maximum node to include.
  --haplotype_order_file HAPLOTYPE_ORDER_FILE
                        File with haplotype order, one per line.
  --haplotype_colors_file HAPLOTYPE_COLORS_FILE
                        TSV file with haplotype and color code.
  --heatmap_order_file HEATMAP_ORDER_FILE
                        File with heatmap sample order, one per line.
```
## Usage - Sample similarity figure
The sample similarity figure can be created by running the chromosome_similarity.py script:
```
python "/mnt/WORKSPACE/psk_workspace/scripts/pangenome/chromosome_similarity.py" -h
usage: chromosome_similarity.py [-h] --input_table INPUT_TABLE --output_dir OUTPUT_DIR --chrom
                                CHROM [--window_size WINDOW_SIZE] --haplotype_colors_file
                                HAPLOTYPE_COLORS_FILE --sample_order_file SAMPLE_ORDER_FILE

Generate chromosome similarity heatmap.

options:
  -h, --help            show this help message and exit
  --input_table INPUT_TABLE
                        Path to the input TSV file with node information.
  --output_dir OUTPUT_DIR
                        Directory to save output files.
  --chrom CHROM         Chromosome to analyze.
  --window_size WINDOW_SIZE
                        Window size for aggregation (default: 200000).
  --haplotype_colors_file HAPLOTYPE_COLORS_FILE
                        TSV file with haplotype and hex color code.
  --sample_order_file SAMPLE_ORDER_FILE
                        File with sample drawing order, one sample per line.
```
## Cite
