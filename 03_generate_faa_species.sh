#!/bin/bash
# Lore Van Santvliet 20/03/2022
# Species tree construction: This script is meant to generate faa files (FASTA amino acid files) from ffn files (FASTA nucleotide of gene regions) for each of the species. These faa files can later be used to generate species trees per species. 

while IFS="," read -r genome species
do
	mkdir -p ../results/intermediate/filtered_faas/"$species"
	transeq ../results/intermediate/filtered_ffns/"$genome".ffn ../results/intermediate/filtered_faas/"$species"/"$genome".faa
	echo $species $genome done
done < <(tail -n +2 ../results/intermediate/genomes_species.csv)
