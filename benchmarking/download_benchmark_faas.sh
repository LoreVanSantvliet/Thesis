#!/bin/bash
# Lore Van Santvliet 07/03/2022
# This script uses entrez direct to download translated CDS GenBank files for the genomes specified in the genome_species_benchmark.csv or genome_species_train.csv file, depending on whether the input argument given is "train" or "benchmark"..

usage=$1

if [ $usage = "training" ]
then
	genomes_species_file=genomes_species_train.csv
else
	genomes_species_file=genomes_species_benchmark.csv
fi	

tail -n+2 ../../results/intermediate/benchmarking/$genomes_species_file | cut -d ',' -f 1 | while read genome
do
	echo ${genome} started
	link=$(esearch -db assembly -query ${genome} < /dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_translated_cds.faa.gz"}')
	wget ${link} -O ${genome}.faa.gz && gunzip -c ${genome}.faa.gz >../../results/intermediate/usage/${genome}"_translated_cds.faa"
	echo ${genome} finished
done

echo reached end
rm *.gz  
