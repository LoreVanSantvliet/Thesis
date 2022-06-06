#!/bin/bash
# Lore Van Santvliet 11/04/2022
# This script removes uninformative positions in the gene tree msas, i.e. positions where many sequences have gaps, using the tool trimAl. 

# script parameters
trim_threshold = 0.1
input_path = ../../results/intermediate/msas
output_path = ../../results/intermediate/msas_trimmed

mkdir $output_path
ls $input_path | while read species
do
	echo "$species"
	if ! test -d $output_path/"${species}"
	then
		echo "$species" " trimmed"
		mkdir $output_path/"$species"
		ls $input_path/"$species" | while read msa_file
		do
			echo $msa_file
			trimal -in $input_path/"${species}"/${msa_file} -out $output_path/"${species}"/${msa_file} -gt $trim_threshold
		done
	fi
	echo "$species" " done"
done
