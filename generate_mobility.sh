#!/bin/bash
# Lore Van Santvliet
# This script is meant to generate mobility files for all genomes.

cat ../results/intermediate/genomes.csv | while read genome
do
	python 02_visualize_mobility_score.py $genome
	echo $genome
done
