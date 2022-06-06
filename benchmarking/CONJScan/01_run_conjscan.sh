#!/bin/bash
# Lore Van Santvliet 07/04/2022
# This script runs ConjScan on multiple genomes.

# method can either be "training", "benchmarking" or "other"
method=$1

input_path=../../../results/intermediate/$method
output_path=../../../results/intermediate/$method/conjscan_raw

mkdir $output_path

mkdir $output_path

for seq in $input_path/*faa; do
	name=${seq##*/}
	name=${name%_translated_cds.faa}
	echo ${name}
        macsyfinder --models conj_chromosome T4SS_typeFATA T4SS_typeFA \
		--db-type ordered_replicon \
		--replicon-topology linear \
                --models-dir Conjugation \
                --sequence-db "$seq"
	mv mac*/all_systems.tsv $output_path/${name}.csv
	rm -r mac*
done

