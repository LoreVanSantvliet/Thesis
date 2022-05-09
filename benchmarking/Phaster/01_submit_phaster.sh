#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script submits jobs to phaster, to detect prophages in 50 randomly selected genomes.

input_path=../../../results/intermediate/benchmarking
output_path=../../../results/intermediate/benchmarking/phaster_raw

mkdir $output_path

tail -n+2 $input_path/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	wget --post-file=$input_path"/"${genome}"_genomic.fna" "http://phaster.ca/phaster_api?contigs=1" -O $output_path/${genome}_job.txt
done