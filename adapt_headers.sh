#!/bin/bash
# Lore Van Santvliet 10/05/2022

# This script is meant to adapt file headers for FASTA files. 

input_path=../results/intermediate/sel_orthogroups/fastas
output_path=../results/intermediate/sel_orthogroups/fastas_adapted
pangenome_path=../results/intermediate/sel_pan/pangenome.tsv

mkdir $output_path

declare -a arr
cat $pangenome_path | while read gene genome orthogroup
do
	echo "$gene"
	echo $genome
	gene="${gene/./_}"
	arr["$gene"]="$genome"
	#echo ${arr["$gene"]}
done

echo ${arr["$gene"]}
echo "start reading fastas"

ls $input_path | while read fasta
do
	echo $fasta
	# touch $output_path/$fasta
	content=$(cat $input_path/$fasta)
	content="${content/./_}"
	#echo $content
	for i in "${!arr[@]}"
	do
		echo $i
		echo arr["$i"]
		content="${content/"$i"/"${arr["$i"]}"}"
	done
	echo $content > $output_path/$fasta
	#grep -o '^[^ ]*' $content | grep .  > $output_path/$fasta
	
done

ls $output_path | while read fasta
do
	grep -o '^[^ ]*' $output_path/$fasta | grep .  > $output_path/$fasta
done
