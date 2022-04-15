#!/bin/bash
# Lore Van Santvliet 04/03/2022
# This script selects a number of genomes out of the genomes_species.csv file, to serve as benchmark genomes for the MGE detection tool. 

n_genomes=50

mkdir ../../results/intermediate/benchmarking
tail -n+2 ../../results/intermediate/genomes_species.csv | gshuf -n $n_genomes > ../../results/intermediate/benchmarking/genomes_species_benchmark.csv
