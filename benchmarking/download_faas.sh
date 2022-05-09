#!/bin/bash
# Lore Van Santvliet 07/03/2022
# This script uses entrez direct to download translated CDS GenBank files for the genomes specified in the genome_species_benchmark.csv file.

tail -n+2 ../../results/intermediate/benchmarking/genomes_species_benchmark.csv | cut -d ',' -f 1 | while read genome
do
	echo ${genome} started
	link=$(esearch -db assembly -query ${genome} < /dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_translated_cds.faa.gz"}')
	wget ${link} -O ${genome}.faa.gz && gunzip -c ${genome}.faa.gz >../../results/intermediate/benchmarking/${genome}"_translated_cds.faa"
	echo ${genome} finished
done

echo reached end
rm *.gz  
