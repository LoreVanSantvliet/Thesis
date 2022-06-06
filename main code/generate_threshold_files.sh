#!/bin/bash
# Lore Van Santvliet 02/05/2022
# This script generates threshold files and sets up performance files for the training of the MGE detection algorithm.

# script parameters
output_path=../../results/intermediate/training

# baseline
echo "10
20
30" > $output_path/baseline_thresholds.txt

# simple
echo -n "" > $output_path/simple_thresholds.txt
for score_threshold in {0.2,0.3,0.4}
do
	for n_threshold in 10 20 30
	do
		echo $score_threshold $n_threshold >> $output_path/simple_thresholds.txt
	done
done

# intermediate
echo -n "" > $output_path/intermediate_thresholds.txt
for score_threshold in {0.2,0.3,0.4}
do
    for n_threshold in 10 20 30
    do
        for initial_n_threshold in 4 5 6
        do
            for distance_threshold in 3 4 5
            do
                echo $score_threshold $n_threshold $initial_n_threshold $distance_threshold >> $output_path/intermediate_thresholds.txt
            done
        done
    done
done


# sophisticated
echo -n "" > $output_path/sophisticated_thresholds.txt
for score_threshold in {0.2,0.3,0.4}
do
    for n_threshold in 10 20 30
    do
        for fraction_threshold in {0.1,0.2,0.3}
        do
            for n_scc_threshold in 3 4 5
            do
                echo $score_threshold $n_threshold $fraction_threshold $n_scc_threshold >> $output_path/sophisticated_thresholds.txt
            done
        done
    done
done
