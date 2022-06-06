#!/bin/bash
# Lore Van Santvliet
# This script is meant to generate mobility files for all genomes.

output_path_files=../../results/mobility_files
output_path_plots=../../results/mobility_plots
output_dir=../../results/intermediate/benchmarking

mkdir $output_path_files/$output_dir
mkdir $output_path_plots/$output_dir

#cat ../results/intermediate/genomes.csv | while read genome
tail -n+2 ../results/intermediate/$output_dir/genomes_species_benchmark.csv | while read genome_species
do
    genome=(${genome_species//,/ })
    echo $genome
	python ../main\code/02_record_genomewide_mobility_scores.py $genome $output_dir
done
