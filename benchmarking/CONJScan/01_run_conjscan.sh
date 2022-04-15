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

#	for conj_type in T4SS_typeFATA T4SS_typeFA
#	do
#        	macsyfinder --models conj_chromosome "$conj_type" \
#                    --db-type ordered_replicon \
#		    --replicon-topology linear \
#                    --models-dir Conjugation \
#                    --sequence-db "$seq"
#		mv mac*/all_systems.tsv ../../results/intermediate/benchmarking/conjscan_raw/${name}-${conj_type}-chrom.csv
#		rm -r mac*
#	done
done

# failed:
#macsyfinder --models conj_plasmid T4SS_typeFATA T4SS_typeFA --db-type ordered_replicon --replicon-topology linear --models-dir Conjugation --sequence-db ../../../results/intermediate/benchmarking/failed/GCA_002092815.1_translated_cds.faa
