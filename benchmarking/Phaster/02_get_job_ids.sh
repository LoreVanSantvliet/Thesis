#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script looks at job submission files of phaster (GCA_..._job.txt files), extracts the job_ids, and pastes them into a new file.

input_path=../../../results/intermediate/benchmarking
output_path=../../../results/intermediate/benchmarking/phaster_raw/job_ids.txt

tail -n+2 $input_path/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	echo $genome
	echo $genome,$(head $input_path/phaster_raw/${genome}_job.txt | cut -d '"' -f 4) >> $output_path
done
