#!/bin/bash
# Lore Van Sanntvliet 20/03/2022
# Species tree construction: This script uses scarap to generate supermatrices for all species, using the .faa files for these species. The supermatrices can later be used to generate species trees. 

# script parameters
input_path = ../results/intermediate/filtered_faas
output_path = ../results/intermediate/supermatrix

ls $input_path | while read species
do
	ls $input_path/"$species"/*.faa > $output_path/faapaths.txt
	scarap core-pipeline $output_path/faapaths.txt core -t 4
	scarap supermatrix $output_path/faapaths.txt core/coregenome.tsv supermatrix
	cp supermatrix/supermatrix_aas.fasta $output_path/"$species".fasta
	rm -r supermatrix
	rm -r core
	rm -r core.prev0
    rm -r core.prev1
	echo $species done
done
rm $output_path/faapaths.txt
