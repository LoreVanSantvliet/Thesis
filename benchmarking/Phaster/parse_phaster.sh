#!/bin/bash
# Lore Van Santvliet 18/04/2022
# This script runs the python script to parse phaster files.

# script parameters
dir=$1
project_path=../../../
input_path=$project_path/results/intermediate/$dir/phaster_raw
output_path=$project_path/results/intermediate/$dir/phaster_parsed

mkdir $output_path

ls $input_path  | while read phaster_file
do
    python ParsePhaster.py $phaster_file $dir
done


