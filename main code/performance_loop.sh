#!/bin/bash

# Lore Van Santvliet 02/05/2022
# Runs the MGE_detection.py and measure_performance.py scripts for a range of thresholds for all training genomes.

# script parameters
directory="benchmarking"
mkdir ../../results/intermediate/$directory
baseline_threshold_path=../../results/intermediate/$directory/baseline_thresholds.txt
simple_threshold_path=../../results/intermediate/$directory/simple_thresholds.txt
intermediate_threshold_path=../../results/intermediate/$directory/intermediate_thresholds.txt
sophisticated_threshold_path=../../results/intermediate/$directory/sophisticated_thresholds.txt
input_path_conjscan=../../results/intermediate/$directory/conjscan_parsed
input_path_phaster=../../results/intermediate/$directory/phaster_parsed
input_path_self_baseline=../../results/MGE_files_baseline/$directory
input_path_self_simple=../../results/MGE_files_simple/$directory
input_path_self_sophisticated=../../results/MGE_files_sophisticated/$directory
input_path_self_intermediate=../../results/MGE_files_intermediate/$directory

method=$1

mkdir ../../results/intermediate/$directory/$method

#rm $file_path_conjscan
#rm $file_path_phaster
#rm $file_path_joined

mob_files_path=../../results/mobility_files/$directory
#mge_files_contig_path=../../results/MGE_files_contig/$directory
#mge_files_path=../../results/MGE_files/$directory
mge_files_contig_baseline_path=../../results/MGE_files_contig_baseline/$directory
mge_files_baseline_path=../../results/MGE_files_baseline/$directory
mge_files_contig_simple_path=../../results/MGE_files_contig_simple/$directory
mge_files_simple_path=../../results/MGE_files_simple/$directory
mge_files_contig_intermediate_path=../../results/MGE_files_contig_intermediate/$directory
mge_files_intermediate_path=../../results/MGE_files_intermediate/$directory
mge_files_contig_sophisticated_path=../../results/MGE_files_contig_sophisticated/$directory
mge_files_sophisticated_path=../../results/MGE_files_sophisticated/$directory

if [ $directory = "training" ]
then
	genomes_file=genomes_train.csv
else
	genomes_file=genomes_benchmark.csv
fi

echo $method

if [ $method = "simple" ]
then
    echo "Simple"
    threshold_path=$simple_threshold_path
    echo $threshold_path
    mge_files_contig_path=$mge_files_contig_simple_path
    mge_files_path=$mge_files_simple_path
    input_path_self=$input_path_self_simple
    
    #echo "genome,n_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_conjscan
    #echo "genome,n_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_phaster
    #echo "genome,n_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_joined
else
    if [ $method = "intermediate" ]
    then
        echo "Intermediate"
        threshold_path=$intermediate_threshold_path
        mge_files_contig_path=$mge_files_contig_intermediate_path
        mge_files_path=$mge_files_intermediate_path
        input_path_self=$input_path_self_intermediate
        
        #echo "genome,score_threshold,n_threshold,initial_n_threshold,distance_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_conjscan
        #echo "genome,score_threshold,n_threshold,initial_n_threshold,distance_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_phaster
        #echo "genome,score_threshold,n_threshold,initial_n_threshold,distance_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_joined    
    else
	if [ $method = "sophisticated" ]
	then        
		echo "Sophisticated"
        	threshold_path=$sophisticated_threshold_path
        	mge_files_contig_path=$mge_files_contig_sophisticated_path
        	mge_files_path=$mge_files_sophisticated_path
        	input_path_self=$input_path_self_sophisticated
        
        	#echo "genome,score_threshold,n_threshold,direction_threshold,n_scc_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_conjscan
        	#echo "genome,score_threshold,n_threshold,direction_threshold,n_scc_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_phaster
        	#echo "genome,score_threshold,n_threshold,direction_threshold,n_scc_threshold,coverage,dispersion,recall,precision,F-score" > $file_path_joined

        else # baseline
		echo "Baseline"
		threshold_path=$baseline_threshold_path
                mge_files_contig_path=$mge_files_contig_baseline_path
                mge_files_path=$mge_files_baseline_path
                input_path_self=$input_path_self_baseline
		
		#echo "genome,score_threshold,n_threshold,coverage,dispersion,recall,precision,F-score">$file_path_conjscan
        	#echo "genome,score_threshold,n_threshold,coverage,dispersion,recall,precision,F-score">$file_path_phaster
        	#echo "genome,score_threshold,n_threshold,coverage,dispersion,recall,precision,F-score">$file_path_joined
	fi
    fi
fi

parameter_code=0
cat $threshold_path | while read threshold_set
do
	echo "Thresholds: " "${threshold_set}"

	#rm -r $mge_files_contig_path
	#rm -r $mge_files_path
	mkdir -p $mge_files_contig_path
	mkdir -p $mge_files_path

	echo "Start detection MGEs"

    	ls $mob_files_path | while read genome_contig_file
    	do
        	python MGE_detection.py $genome_contig_file $method "$threshold_set" "$threshold_set"
        	echo $genome_contig_file " done"
    	done
    
    	echo "Start merging contig MGE files"
   	rm -r $mge_files_path/"$threshold_set" 
    	#mkdir -p $mge_files_path
    
	ls $mge_files_contig_path/"$threshold_set" | while read genome_contig_file
    	do
            	genome=${genome_contig_file%-[A-Z]*[0-9]*.[0-9]*.csv}
            	genome_file=$genome".csv"
		mkdir -p $mge_files_path/"$threshold_set"
            	if [[ ! -e $mge_files_path/"$threshold_set"/$genome_file ]]
            	then
                    	cp $mge_files_contig_path/"$threshold_set"/$genome_contig_file $mge_files_path/"$threshold_set"/$genome_file
            	else
                    	tail -n +2 $mge_files_contig_path/"$threshold_set"/$genome_contig_file >> $mge_files_path/"$threshold_set"/$genome_file
            	fi
            	echo $genome_file " merged"
    	done

    	echo "Start performance measurements"
  	echo $input_path_self 
	
	file_path_conjscan=../../results/intermediate/$directory/$method/"$threshold_set"/conjscan_performance.csv
	file_path_phaster=../../results/intermediate/$directory/$method/"$threshold_set"/phaster_performance.csv
	file_path_joined=../../results/intermediate/$directory/$method/"$threshold_set"/joined_performance.csv
 	mkdir ../../results/intermediate/$directory/$method/"$threshold_set"
	rm "$file_path_conjscan"
	echo "conjscan file removed"
	rm "$file_path_phaster"
	echo "phaster file removed"
	rm "$file_path_joined"
	echo "joined file removed"
	echo "$threshold_code"
    	cat ../results/intermediate/$directory/$genomes_file | while read genome_file
	do
        	echo "start comparison " $genome_file
        	python measure_performance.py $method $input_path_self/"$threshold_set"/$genome_file $input_path_conjscan/$genome_file "$file_path_conjscan" $input_path_phaster/$genome_file "$file_path_phaster" "$file_path_joined" $(basename $genome_file .csv) "$threshold_set"
        	echo "comparison " $genome_file " finished"
    	done
done


