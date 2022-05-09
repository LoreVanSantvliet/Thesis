#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script selects a number of genomes out of the genomes_species.csv file, to serve as parameter finetuning and benchmark genomes for the MGE detection tool. 

n_genomes=50
input_path=../../results/intermediate/genomes_species.csv
benchmarking_path=../../results/intermediate/benchmarking
training_path=../../results/intermediate/training

mkdir $benchmarking_path
mkdir $training_path
tail -n+2 $input_path | gshuf -n $((n_genomes*2)) > $benchmarking_path/genomes_species_benchmark_int.csv

echo "genome,species" > $training_path/genomes_species_train.csv
echo "genome,species" > $benchmarking_path/genomes_species_benchmark.csv

# first 50 files --> training
head -n $n_genomes $benchmarking_path/genomes_species_benchmark_int.csv >> $training_path/genomes_species_train.csv
# last 50 files --> benchmarking
tail -n+$(($n_genomes+1)) $benchmarking_path/genomes_species_benchmark_int.csv >> $benchmarking_path/genomes_species_benchmark.csv

rm $benchmarking_path/genomes_species_benchmark_int.csv