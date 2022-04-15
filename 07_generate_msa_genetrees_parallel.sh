#!/bin/bash
# Lore Van Santvliet 30/03/2022
# This script generates multiple sequence alignments (MSAs) for each orthogroup in the ../results/intermediate/gene_tree_sequences directory, using MAFFT, based on the orthogroup FASTA files, per species, in a parallel fashion.

# script parameters
n_threads=12

make_aln(){
	#echo $2
	orthogroup=$1
	species="$2"
	species_2=$3
	n_seq=$(grep '>' ../results/intermediate/gene_tree_sequences/${species}\ ${species_2}/$orthogroup | wc -l)
        if [[ $n_seq -gt 4 ]]
        then	
                echo $orthogroup
                mafft  --quiet --thread 8 ../results/intermediate/gene_tree_sequences/${species}\ ${species_2}/$orthogroup >  ../results/intermediate/msas/${species}\ ${species_2}/$orthogroup
        mafft  --quiet --thread 8 ../results/intermediate/gene_tree_sequences/${species}\ ${species_2}/"$orthogroup" >  ../results/intermediate/msas/${species}\ ${species_2}/"$orthogroup"
	fi
}

export -f make_aln

#ls ../results/intermediate/gene_tree_sequences | while read species
less ../results/intermediate/progress | while read species
do
	echo "$species"
    	mkdir ../results/intermediate/msas/"$species"
    ls ../results/intermediate/gene_tree_sequences/"$species" | parallel \
    	--jobs $n_threads \
	--no-notice \
	--verbose \
	make_aln {} "$species"
    echo $species " done"
done
