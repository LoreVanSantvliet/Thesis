#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script uses entrez direct to download genomic GenBank files for the genomes specified in the genome_species_benchmark.csv file.

tail -n+2 ../../results/intermediate/benchmarking/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	echo ${genome} started
	link=$(esearch -db assembly -query ${genome} < /dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}')
	wget ${link} -O ${genome}.fna.gz && gunzip -c ${genome}.fna.gz >../../results/intermediate/benchmarking/${genome}"_genomic.fna"
	echo ${genome} finished
done

echo reached end
rm *.gz  
