#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script looks at job submission files of phaster (GCA_..._job.txt files), extracts the job_ids, and pastes them into a new file.

tail -n+2 ../../../results/intermediate/benchmarking/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	echo $genome
	echo $genome,$(head ../../../results/intermediate/benchmarking/phaster/${genome}_job.txt | cut -d '"' -f 4) >> ../../../results/intermediate/benchmarking/phaster_raw/job_ids.txt
done
