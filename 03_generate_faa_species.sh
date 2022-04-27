#!/bin/bash
# Lore Van Santvliet 20/03/2022
# Species tree construction: This script is meant to generate faa files (FASTA amino acid files) from ffn files (FASTA nucleotide of gene regions) for each of the species. These faa files can later be used to generate species trees per species. 

# script parameters
input_path = ../results/intermediate/filtered_ffns
output_path = ../results/intermediate/filtered_faas
genomes_species_file = ../results/intermediate/genomes_species.csv

while IFS="," read -r genome species
do
	mkdir -p $output_path/"$species"
	transeq $input_path/"$genome".ffn $output_path/"$species"/"$genome".faa
	echo $species $genome done
done < <(tail -n +2 genomes_species_file)
