#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script submits jobs to phaster, to detect prophages in 50 randomly selected genomes.

tail -n+2 ../../../results/intermediate/benchmarking/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	wget --post-file="../../../results/intermediate/benchmarking/"${genome}"_genomic.fna" "http://phaster.ca/phaster_api?contigs=1" -O ../../../results/intermediate/benchmarking/phaster/${genome}_job.txt	
done
