#!/bin/bash
# Lore Van Santvliet 11/04/2022
# This script changes tree node labels from gene_contig format to gene format.

ls ../results/intermediate/gene_trees | while read species
do
	sed 's/-[A-Z]*[0-9]*\.[0-9]*_[0-9]*//g' ../results/intermediate/"$species"/*.nw > ../results/intermediate/gene_trees_mod/"$species"/*.nw
done
