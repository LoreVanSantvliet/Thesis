#!/bin/bash
# Lore Van Santvliet 31/03/2022
# This script removes '\\n' characters in a file (e.g. Phaster output file) and replaces them by '\n', a newline character. It also replaces a line with '-' characters by spaces.

input_path=../../../results/intermediate/benchmarking/phaster_raw

ls $input_path | while read phaster_file
do
	sed -i '' 's/\\n/\'$'\n/g' $input_path/$phaster_file
	sed -i '' 's/--/  /g' $input_path/$phaster_file
done

