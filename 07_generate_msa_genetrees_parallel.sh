#!/bin/bash
# Lore Van Santvliet 30/03/2022
# This script generates multiple sequence alignments (MSAs) for each orthogroup in the ../results/intermediate/gene_tree_sequences directory, using MAFFT, based on the orthogroup FASTA files, per species, in a parallel fashion.

# script parameters
n_threads=12
input_path = ../results/intermediate/gene_tree_sequences
output_path = ../results/intermediate/msas

make_aln(){
	orthogroup=$1
	species="$2"
	species_2=$3
	n_seq=$(grep '>' $input_path/${species}\ ${species_2}/$orthogroup | wc -l)
        if [[ $n_seq -gt 4 ]]
        then	
                echo $orthogroup
                mafft  --quiet $input_path/${species}\ ${species_2}/$orthogroup >  $output_path/${species}\ ${species_2}/$orthogroup
        mafft  --quiet $input_path/${species}\ ${species_2}/"$orthogroup" >  $output_path/${species}\ ${species_2}/"$orthogroup"
        else
            echo $orthogroup >> ../results/intermediate/no_gene_tree.csv
            
	fi
}

export -f make_aln

less .$input_path | while read species
do
	echo "$species"
    	mkdir $output_path/"$species"
    ls input_path/"$species" | parallel \
    	--jobs $n_threads \
	--no-notice \
	--verbose \
	make_aln {} "$species"
    echo $species " done"
done
