#!/bin/bash
# Lore Van Santvliet 07/04/2022
# This script runs ConjScan on multiple genomes. It needs 

for seq in ../../../results/intermediate/benchmarking/*faa; do
	name=${seq##*/}
	name=${name%_translated_cds.faa}
	echo ${name}
        macsyfinder --models conj_chromosome T4SS_typeFATA T4SS_typeFA \
		--db-type ordered_replicon \
		--replicon-topology linear \
                --models-dir Conjugation \
                --sequence-db "$seq"
	mv mac*/all_systems.tsv ../../../results/intermediate/benchmarking/conjscan_raw/${name}.csv
	rm -r mac*
done

