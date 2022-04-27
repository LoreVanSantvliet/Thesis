#!/bin/bash
# Lore Van Santvliet 16/04/2022
# This script automatically runs the MGE_detection.py script for the available mobility frames. It also pastes the different files belonging to the same genome together.

#ls ../results/mobility_files | while read genome_contig_file
#do
#	python MGE_detection.py $genome_contig_file
#	echo $genome_contig_file " done"
#done

mkdir ../results/MGE_files

ls ../results/MGE_files_contig | while read genome_contig_file
do
	genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
	genome_file=$genome".csv"
	if [[ ! -e ../results/MGE_files/$genome_file ]]
	then
		cp ../results/MGE_files_contig/$genome_contig_file ../results/MGE_files/$genome_file
	else
		tail -n +2 ../results/MGE_files_contig/$genome_contig_file >> ../results/MGE_files/$genome_file
	fi
	echo $genome_file
done


mkdir ../results/MGE_files_simple

ls ../results/MGE_files_contig_simple | while read genome_contig_file
do
        genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
        genome_file=$genome".csv"
        if [[ ! -e ../results/MGE_files/$genome_file ]]
        then
                cp ../results/MGE_files_contig/$genome_contig_file ../results/MGE_files_simple/$genome_file
        else
                tail -n +2 ../results/MGE_files_contig_simple/$genome_contig_file >> ../results/MGE_files_simple/$genome_file
        fi
        echo $genome_file
done

