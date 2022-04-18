#!/bin/bash
# Lore Van Santvliet 31/03/2022
# This script removes '\\n' characters in a file (e.g. Phaster output file) and replaces them by '\n', a newline character. It also replaces a line with '-' characters by spaces.

ls ../../../results/intermediate/benchmarking/phaster | while read phaster_file
do
	sed -i '' 's/\\n/\'$'\n/g' ../../../results/intermediate/benchmarking/phaster/$phaster_file
	sed -i '' 's/--/  /g' ../../../results/intermediate/benchmarking/phaster/$phaster_file
done

