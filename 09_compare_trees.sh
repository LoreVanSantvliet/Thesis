#!/bin/bash
# Lore Van Santvliet 22/03/2022
# This script is meant to compare gene trees and species trees using ete3 (treeKO). 

# script parameters
n_threads=8

species="Apilactobacillus kunkeei_A"
ete3 compare -cpu $n_threads -t ../results/intermediate/gene_trees/"$species"/F00007_3.txt -r ../results/intermediate/species_trees/"$species".fasta --unrooted --treeko


#ete3 compare -cpu $n_threads -t ../results/intermediate/gene_trees/"$species"/F00007_3.txt ../results/intermediate/gene_trees/"$species"/F00008_2.txt -r ../data/species_trees/"$species".fasta --unrooted --treeko

#mkdir ../data/comparisons/"$species"
#ete3 compare -cpu $n_threads -t ../results/intermediate/gene_trees/"$species"/F00007_3.txt ../results/intermediate/gene_trees/"$species"/F00008_2.txt -r ../results/intermediate/species_trees/"$species".fasta --unrooted --treeko > ../results/intermediate/comparisons/"$species"/F00008_2.txt 


#ls ../data/gene_trees | while read species
#do
	#mkdir ../results/intermediate/comparisons/"$species"
	#ls ../results/intermediate/gene_trees/"$species" | while read orthogroup_file
	#do
		#ete3 compare -cpu $n_threads -t ../results/intermediate/gene_trees/"$species"/"orthogroup_file" -r ../results/intermediate/species_trees/"$species".fasta --unrooted --treeko > ../results/intermediate/comparisons/"$species"/"$orthogroup_file" 
