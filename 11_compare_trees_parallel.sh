#!/bin/bash
# Lore Van Santvliet 20/04/2022
# This script compares gene trees for all orthogroups and species, to the species trees of those species, in a parallel fashion. It uses the compare_tree.py script for this.

# script parameters
path_gene_trees=../results/intermediate/gene_trees_mod
output_path=../results/intermediate/tree_mobility_frames
n_threads=8

mkdir $output_path
ls $path_gene_trees | while read species
do
	echo "$species"
    echo "orthogroup,treeko_dist" > $output_path/"$species".csv
    ls $path_gene_trees/"$species" | parallel \
    	--jobs $n_threads \
	--no-notice \
	--verbose \
	python compare_trees.py "$species"
    echo "$species" " done"
done
