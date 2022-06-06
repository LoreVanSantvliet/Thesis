#!/bin/bash
# Lore Van Santvliet 12/05/2022
# This script generates threshold files and sets up performance files for the benchmarking of the MGE detection algorithm.

# script parameters
output_path=../../results/intermediate/benchmarking

mkdir $output_path

# simple
echo "10" > $output_path/simple_thresholds.txt

# intermediate
echo "0.4 30 5 4" > $output_path/intermediate_thresholds.txt


# sophisticated
echo "0.4 30 0.2 5" > $output_path/sophisticated_thresholds.txt

