#!/bin/bash
# Lore Van Santvliet 04/03/2022

less ../../../results/intermediate/benchmarking/phaster_raw/job_ids.txt | while read genome_job_id
do
	genome=$(echo ${genome_job_id} | cut -d ',' -f 1)
	job_id=$(echo ${genome_job_id} | cut -d ',' -f 2)
	wget "http://phaster.ca/phaster_api?acc=${job_id}" -O ../../../results/intermediate/benchmarking/phaster_raw/${genome}.txt
done


