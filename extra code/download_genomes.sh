#!/bin/bash
# Lore Van Santvliet 15/03/2022
# This script is meant for downloading genome files from the UA server onto my personal computer.

cat ../../results/intermediate/genomes.csv | while read genome
do
	scp lorevs@143.129.141.185:/media/ssdsata/stijn/projects_phd/legen/results/all/genes/ffns/"$genome".ffn.gz ../../data/filtered_ffns/
	echo $genome
done
