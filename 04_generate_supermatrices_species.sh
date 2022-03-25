#!/bin/bash
# Lore Van Sanntvliet 20/03/2022
# Species tree construction: This script uses scarap to generate supermatrices for all species, using the faa.files for these species. The supermatrices can be used to generate species trees. 

ls ../results/intermediate/filtered_faas | while read species
do
	ls ../results/intermediate/filtered_faas/"$species"/*.faa > ../results/intermediate/supermatrix/faapaths.txt
	$scarap core-pipeline ../results/intermediate/supermatrix/faapaths.txt core -t 4
	$scarap supermatrix ../results/intermediate/supermatrix/faapaths.txt core/coregenome.tsv supermatrix
	cp supermatrix/supermatrix_aas.fasta ../results/intermediate/supermatrix/"$species".fasta
	rm -r supermatrix
	rm -r core
	rm -r core.prev0
    rm -r core.prev1
	echo $species done
done
rm ../results/intermediate/supermatrix/faapaths.txt
