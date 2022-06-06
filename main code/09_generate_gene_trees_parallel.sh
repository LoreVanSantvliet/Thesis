#!/bin/bash
# Lore Van Santvliet 06/04/2022
# This script creates gene trees for all orthogroups in all species. It uses a function in the python script generate_genetrees.py, and speeds up the process by parallellization.

# script parameters
n_threads = 8
input_path = ../../results/intermediate/msas_sel

make_tree(){
	orthogroup=$1
	species=$2
	species_2=$3
	echo "species: " $species $species_2 "    - orthogroup: " $orthogroup
	python3 generate_gene_trees.py $orthogroup $species $species_2
}

export -f make_tree

ls $input_path | while read species
do
	ls $input_path/"$species" | parallel \
		--jobs $n_threads \
		--no-notice \
		--verbose \
		make_tree {} "$species"
done
