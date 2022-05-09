#!/bin/bash
# Lore Van Santvliet 16/04/2022
# This script automatically runs the MGE_detection.py script for the available mobility frames. It also pastes the different files belonging to the same genome together.

# script parameters
mob_files_path=../results/mobility_files/training
mge_files_contig_path=../results/MGE_files_contig/training
mge_files_path=../results/MGE_files/training
mge_files_contig_simple_path=../results/MGE_files_contig_simple/training
mge_files_simple_path=../results/MGE_files_simple/training
mge_files_contig_intermediate_path=../results/MGE_files_contig_intermediate/training
mge_files_intermediate_path=../results/MGE_files_intermediate/training

ls $mob_files_path | while read genome_contig_file
do
	python MGE_detection.py $genome_contig_file "simple" 20
	echo $genome_contig_file " done"
done

#python MGE_detection.py "GCA_000192185.1-GL872376.1.csv" "simple" 20

mkdir $mge_files_path

ls $mge_files_contig_path | while read genome_contig_file
do
    genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
    genome_file=$genome".csv"
    if [[ ! -e $mge_files_path/$genome_file ]]
    then
        cp $mge_files_contig_path/$genome_contig_file $mge_files_path/$genome_file
    else
        tail -n +2 $mge_files_contig_path/$genome_contig_file >> $mge_files_path/$genome_file
    fi
    echo $genome_file
done


mkdir $mge_files_simple_path

ls $mge_files_contig_simple_path | while read genome_contig_file
do
        genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
        genome_file=$genome".csv"
        if [[ ! -e $mge_files_simple_path/$genome_file ]]
        then
                cp $mge_files_contig_simple_path/$genome_contig_file $mge_files_simple_path/$genome_file
        else
                tail -n +2 $mge_files_contig_simple_path/$genome_contig_file >> $mge_files_simple_path/$genome_file
        fi
        echo $genome_file
done

mkdir $mge_files_intermediate_path

ls $mge_files_contig_intermediate_path | while read genome_contig_file
do
        genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
        genome_file=$genome".csv"
        if [[ ! -e $mge_files_intermediate_path/$genome_file ]]
        then
                cp $mge_files_contig_interemdiate_path/$genome_contig_file $mge_files_intermediate_path/$genome_file
        else
                tail -n +2 $mge_files_contig_intermediate_path/$genome_contig_file >> $mge_files_intermediate_path/$genome_file
        fi
        echo $genome_file
done

