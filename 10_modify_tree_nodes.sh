#!/bin/bash
# Lore Van Santvliet 11/04/2022
# This script changes tree node labels: inner leaf names are removed from all gene and species trees, and node names are changed from gene-contig format to gene format for gene trees where this is necessary.

# script parameters
input_path_genes=../results/intermediate/gene_trees
output_path_genes=../results/intermediate/gene_trees_mod
input_path_species=../results/intermediate/species_trees
output_path_species=../results/intermediate/species_trees_mod

mkdir $output_path_genes
mkdir $output_path_species

ls $input_path_genes | while read species
do
	
	mkdir $output_path_genes/"$species"
	ls $input_path_genes/"$species" | while read treefile
	do
		# remove inner node names
		sed 's/Inner[0-9]*//g;s/-[A-Z]*[0-9]*\.[0-9]*_[0-9]*//g' $input_path_genes/"$species"/"$treefile" > $output_path_genes/"$species"/$treefile
	done
	
	sed 's/Inner[0-9]*//g' $input_path_species/"$species".nw > $output_path_species/"$species".nw

done



