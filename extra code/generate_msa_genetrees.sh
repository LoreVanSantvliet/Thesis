#!/bin/bash
# Lore Van Santvliet 16/03/2022
# This script generates multiple sequence alignments (MSAs) for each orthogroup in the ../../results/intermediate/gene_tree_sequences directory, using MAFFT, ased on the orthigroup FASTA files, per species.

ls ../../results/intermediate/gene_tree_sequences | while read species
do
	echo $species " start"
    mkdir ../../results/intermediate/msas/"$species"
    ls ../../results/intermediate/gene_tree_sequences/"$species" | while read orthogroup
    do
        mafft  --quiet --thread 8 ../../results/intermediate/gene_tree_sequences/"$species"/"$orthogroup" >  ../results/intermediate/msas/"$species"/"$orthogroup"
        echo $orthogroup
    done
    echo $species " done"
done
