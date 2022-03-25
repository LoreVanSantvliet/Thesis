#!/bin/bash
#Lore Van Santvliet 23/03/2022
# This script is meant for generating directories corresponding to all species of the analysis.

ls ../results/intermediate/filtered_faas | while read species
do
	mkdir ../results/intermediate/gene_tree_sequences/"$species"
done
