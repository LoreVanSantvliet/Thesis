#!/bin/bash
# Lore Van Santvliet 07/03/2022
# This script uses entrez direct to download translated CDS GenBank files for the genomes specified in the genome_species_benchmark.csv file.

method=$1
echo $method

if [ $method == "benchmarking" ]
then
	genomes_species_file="genomes_species_benchmark.csv"
elif [ $method == "training" ]
then
	genomes_species_file="genomes_species_train.csv"
else
	genomes_species_file="genomes_species_other.csv"
	echo $method
fi

#cat ../../results/intermediate/$method/"$genomes_species_file"

tail -n+2 ../../results/intermediate/$method/"$genomes_species_file" | cut -d ',' -f 1 | while read genome
do
	echo ${genome} started
	link=$(esearch -db assembly -query ${genome} < /dev/null | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_translated_cds.faa.gz"}')
	wget ${link} -O ${genome}.faa.gz && gunzip -c ${genome}.faa.gz >../../results/intermediate/$method/${genome}"_translated_cds.faa"
	echo ${genome} finished
done

echo reached end
rm *.gz  
