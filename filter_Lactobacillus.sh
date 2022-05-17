#!/bin/bash
# Lore Van Santvliet 10/05/2022
# This script filters the genomes belonging to the species Lactobacillus gasseri, from the genomes.metadata.csv file.

project_path=..
input_path=$project_path/data/genomes_metadata.csv
output_path=$project_path/results/intermediate/Lactobacillus_genomes.scv

awk -F "," '$2=="Lactobacillus gasseri" {print $1}' $input_path > $output_path

