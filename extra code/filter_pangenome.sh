#!/bin/bash
# Lore Van Santvliet 20/05/2022

# This script filters the pangenome to only contain genomes from Lactobacillus gasseri.

project_path=../..
pangenome_path=$project_path/data/pangenome.tsv
genome_path=$project_path/results/intermediate/sel_genomes.csv
output_path=$project_path/results/intermediate/sel_pan

scarap filter -g $genome_path $pangenome_path $output_path
