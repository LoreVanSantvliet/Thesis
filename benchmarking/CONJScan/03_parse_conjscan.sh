#!/bin/bash
# Lore Van Santvliet 13/04/2022
# This script removes the first three lines from conjscan output files, if this is required, and then it parses the file using the ParseConjScan.py script.

# script parameters
method=$1 # training, benchmarking or other
input_path=../../../results/intermediate/$method/conjscan_raw
output_path=../../../results/intermediate/$method/conjscan_parsed

mkdir $output_path

ls $input_path | while read conjscan_file
do
	if read -n1 char <$input_path/$conjscan_file; [[ $char = "#" ]]
	then
		sed -i'' -e '1,3d' $input_path/$conjscan_file
		echo "Trimmed " $conjscan_file
	else
		echo $conjscan_file " was already trimmed"
	fi
	python ParseConjScan.py $conjscan_file $method
	echo $conjscan_file " parsed"
done

rm $input_path/*e
echo "Finished"

