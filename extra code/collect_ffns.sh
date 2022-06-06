#!/bin/bash
# Lore Van Santvliet 10/05/2022

# This script collects ffn files and puts them into a directory. (Lactobacilllus gasseri)

project_path=../..
input_path=$project_path/results/intermediate/sel_genomes.csv
ffn_path=$project_path/results/intermediate/filtered_ffns
output_path=$project_path/results/intermediate/sel_ffns

mkdir $output_path

cat $input_path | while read genome
do
	cp $ffn_path/$genome".ffn" $output_path
done

