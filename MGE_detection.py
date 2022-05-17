#!/usr/bin/env python
# coding: utf-8

"""Performs MGE detection in genomes.

Uses the phyletic distribution-based mobility scores of orthogroups to detect MGEs (mobile genetic elements) in contigs. The MGEs are returned in tabular format, where all gene numbers (as appearing in the contig of choice) belonging to a particular MGE are listed. The way this script sould be run is as follows: python MGE_detection.py <genome>-<contig>.csv, where the <genome>-<contig>.csv file is present in the input path (see script parameters).

Input files:
<genome>-<contig>.csv: this is a csv file containing the following columns: order (order of this particular gene in this contig), orthogroup (orthogroup ID), gene (gene ID), accessory_fraction (mobility score), count (how many times does this orthogroup appear in this genome), accessory (is this particular gene accessory in the species to which the genome belongs). 

Output files:
<genome>-<contig>.csv: this is a csv file containing the following columns: contig, MGE (integer indicating the MGE to which these particular gene belongs), gene_nr.
"""

# import statements
import pandas as pd
from pathlib import Path
import numpy as np
import os
import sys

# script parameters
parameter_code=sys.argv[3]
directory="training"
project_path = Path().resolve().parent
input_path = project_path / "results" / "mobility_files" / directory
output_path = project_path / "results" / "MGE_files_contig" / directory / parameter_code
output_path_baseline = project_path / "results" / "MGE_files_contig_baseline" / directory / parameter_code
output_path_simple = project_path / "results" / "MGE_files_contig_simple" / directory / parameter_code
output_path_intermediate = project_path / "results" / "MGE_files_contig_intermediate" / directory / parameter_code
output_path_sophisticated = project_path / "results" / "MGE_files_contig_sophisticated" / directory / parameter_code

# import data
genome_contig_file=sys.argv[1]
print("file: "+genome_contig_file)
frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
method=sys.argv[2]

def end_detection(mob_frame, score_threshold, n_threshold, fraction_threshold, n_scc_threshold):
    # The index of the end boundary gene of the MGE equals the index of the last gene not belonging to the MGE (returned by applying the get_indices_lengths function on the mobility frame in reverse order) plus the length of the MGE. The gene_nr of thi bundary equals this index +1.
    genes,lengths=get_indices_lengths(mob_frame[::-1], score_threshold, n_threshold, fraction_threshold, n_scc_threshold)
    end_genes=list(np.array(genes)+np.array(lengths))
    return [gene+1 for gene in end_genes][::-1], lengths[::-1]

def start_detection(mob_frame, score_threshold, n_threshold, fraction_threshold, n_scc_consecutive):
    # The index of the start boundary gene of the MGE equals the index of the first gene not belonging to the MGE (returned by applying the get_indices_lengths function on the mobility frame) minus the length of the MGE. The gene_nr of this bundary equals this index +1.
    genes,lengths=get_indices_lengths(mob_frame, score_threshold, n_threshold, fraction_threshold, n_scc_consecutive)
    start_genes=list(np.array(genes)-np.array(lengths))
    return [gene+1 for gene in start_genes], lengths

def get_indices_lengths(frame, score_threshold, n_threshold, fraction_threshold, n_scc_threshold):
    """Returns the index of the first gene not belonging to an MGE anymore when looping through consecutive genes in a contig, and length of the MGE."""
    i=0
    exc=0
    n_scc_consecutive=0
    indices=[]
    lengths=[]
    for index,row in frame.iterrows():
        #consecutive=(row.accessory_fraction>score_threshold)
        if row.accessory_fraction>=float(score_threshold):
            # look for consecutive genes with an elevated mobility score
            i=i+1
        elif row.accessory|(row['count']>1):
            #exc are exceptions, so genes that do not have an elevated mobility score
            exc=exc+1
            # in case these exceptions are multi-copy or accessory, or when they do not appear too close to the start of the potential MGE, they are tolerated
            i=i+1
        elif ((exc<float(fraction_threshold)*i)&(n_scc_consecutive<int(n_scc_threshold))):
            exc=exc+1
            # in case these exceptions do not appear too often and not too often in a consecutive fashion, they are also tolerated
            i=i+1
            n_scc_consecutive = n_scc_consecutive + 1
        elif i>=int(n_threshold):
            # if the number of consecutive genes that belong to a potential MGE is large enough, the potential MGE is predicted to be a true MGE
            # in this case, the index and length of the MGE are saved
            indices.append(index)
            lengths.append(i)
            # parameters are set to zero again, to find a new (potential) MGE
            i=0
            exc=0
            n_scc_consecutive = 0
        else:
            # this case captures everything else: genes that are not part of an MGEs
            # parameters are set to zero again, to find a new (potential) MGE
            i=0
            exc=0
            n_scc_consecutive = 0
    return indices,lengths

def detect_mges_baseline(frame, n_threshold):
    i=0
    start_genes=[]
    end_genes=[]
    for index,row in frame.iterrows():
        consecutive=row.accessory|(row['count']>1)
        if consecutive:
            i=i+1
        else:
            if i>=int(n_threshold):
                start_genes.append(index+1-i)
                end_genes.append(index)
                i=0
            else:
                i=0
    return start_genes, end_genes

def detect_mges_simple(frame, score_threshold, n_threshold):
    i=0
    start_genes=[]
    end_genes=[]
    for index,row in frame.iterrows():
        consecutive=row.accessory_fraction>float(score_threshold)
        if consecutive:
            i=i+1
        else:   
            if i>=int(n_threshold):
                start_genes.append(index+1-i)
                end_genes.append(index)
                i=0
            else:
                i=0
    return start_genes, end_genes

def detect_mges_initial(frame, score_threshold, initial_n_threshold):
    i=0
    start_genes=[]
    end_genes=[]
    for index,row in frame.iterrows():
        consecutive=row.accessory_fraction>float(score_threshold)
        if consecutive:
            i=i+1
        else:
            if i>=int(initial_n_threshold):
                start_genes.append(index+1-i)
                end_genes.append(index)
                i=0
            else:
                i=0
    return start_genes, end_genes

def join_mges_initial(start_genes, end_genes, n_threshold, distance_threshold):
    indices_to_join=[]
    for index in range(len(start_genes)):
        if index != 0:
            if(start_genes[index]<=end_genes[index-1]+int(distance_threshold)):
                indices_to_join.append(index)
    
    for index in range(len(indices_to_join))[::-1]:
        start_genes.remove(start_genes[index])
        end_genes[index]=end_genes[index-1]
        end_genes.remove(end_genes[index-1])
        
    for index in range(len(start_genes))[::-1]:
        if end_genes[index]-start_genes[index]<int(n_threshold):
            start_genes.remove(start_genes[index])
            end_genes.remove(end_genes[index])
    
    return start_genes, end_genes

def detect_mges_intermediate(frame, score_threshold, n_threshold, initial_n_threshold, distance_threshold):
    start_genes, end_genes=detect_mges_initial(frame, score_threshold, initial_n_threshold)
    start_genes, end_genes=join_mges_initial(start_genes, end_genes, n_threshold, distance_threshold)
    return start_genes, end_genes
    
def generate_output(genome_contig_file, score_threshold, n_threshold, fraction_threshold, n_scc_threshold):
    """Generates csv files, containing a columns with the contig, a column with a number indicating an individual MGE, and a columns with gene_nrs, indicating the genes that belong to a particular MGE."""
    contig_extension=genome_contig_file.split('-')[1]
    contig=str(contig_extension.split('.')[0]+'.'+contig_extension.split('.')[1])
    MGE_list=[]
    mob_frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
    start_genes, start_lengths=start_detection(mob_frame, score_threshold, n_threshold, fraction_threshold, n_scc_threshold)
    end_genes, end_lengths=end_detection(mob_frame, score_threshold, n_threshold, fraction_threshold, n_scc_threshold)
    
    MGE=0
    if (len(start_genes)>0)&(len(end_genes)>0):
        for i in range(0,len(start_genes)):
            for j in range(0, len(end_genes)):
                if (start_genes[i]+start_lengths[i] >= end_genes[j])&(end_genes[j]-end_lengths[j] <= start_genes[i])&(start_genes[i]<end_genes[j]):
                    MGE=MGE+1
                    for gene in range(start_genes[i], end_genes[j]+1):
                        MGE_list.append({'contig':contig,'MGE':MGE,'gene_nr':gene})
    
    #if (len(start_genes)>0)&(len(end_genes)>0):
    #    if(len(start_genes)!=len(end_genes)):
    #        print("PROBLEM: start and end boundary mismatch")
    #    for i in range(0,len(start_genes)):
    #        for gene in range(start_genes[i], end_genes[i]+1):
    #            MGE_list.append({'contig':contig,'MGE':i+1,'gene_nr':gene})
    
    MGE_frame = pd.DataFrame.from_records(MGE_list)
    if(not MGE_frame.empty):
        MGE_frame.to_csv(output_path_sophisticated / genome_contig_file) 

def generate_output_baseline(genome_contig_file,  n_threshold):
    """Generates csv files, containing a columns with the contig, a column with a number indicating an individual MGE, and a columns with gene_nrs, indicating the genes that belong to a particular MGE."""
    contig_extension=genome_contig_file.split('-')[1]
    contig=str(contig_extension.split('.')[0]+'.'+contig_extension.split('.')[1])
    MGE_list=[]
    mob_frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
    start_genes, end_genes=detect_mges_baseline(mob_frame, n_threshold)
    for i in range(0,len(start_genes)):
        for gene in range(start_genes[i], end_genes[i]+1):
            MGE_list.append({'contig':contig,'MGE':i+1,'gene_nr':gene})
    MGE_frame = pd.DataFrame.from_records(MGE_list)
    if(not MGE_frame.empty):
        MGE_frame.to_csv(output_path_baseline / genome_contig_file) 
        
def generate_output_simple(genome_contig_file, score_threshold, n_threshold):
    """Generates csv files, containing a columns with the contig, a column with a number indicating an individual MGE, and a columns with gene_nrs, indicating the genes that belong to a particular MGE."""
    contig_extension=genome_contig_file.split('-')[1]
    contig=str(contig_extension.split('.')[0]+'.'+contig_extension.split('.')[1])
    MGE_list=[]
    mob_frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
    start_genes, end_genes=detect_mges_simple(mob_frame, score_threshold, n_threshold)
    for i in range(0,len(start_genes)):
        for gene in range(start_genes[i], end_genes[i]+1):
            MGE_list.append({'contig':contig,'MGE':i+1,'gene_nr':gene})
    MGE_frame = pd.DataFrame.from_records(MGE_list)
    if(not MGE_frame.empty):
        MGE_frame.to_csv(output_path_simple / genome_contig_file)


def generate_output_intermediate(genome_contig_file, score_threshold, n_threshold, initial_n_threshold, distance_threshold):
    """Generates csv files, containing a columns with the contig, a column with a number indicating an individual MGE, and a columns with gene_nrs, indicating the genes that belong to a particular MGE."""
    contig_extension=genome_contig_file.split('-')[1]
    contig=str(contig_extension.split('.')[0]+'.'+contig_extension.split('.')[1])
    MGE_list=[]
    mob_frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
    start_genes, end_genes=detect_mges_intermediate(mob_frame, score_threshold, n_threshold, initial_n_threshold, distance_threshold)
    for i in range(0,len(start_genes)):
        for gene in range(start_genes[i], end_genes[i]+1):
            MGE_list.append({'contig':contig,'MGE':i+1,'gene_nr':gene})
    MGE_frame = pd.DataFrame.from_records(MGE_list)
    if(not MGE_frame.empty):
        MGE_frame.to_csv(output_path_intermediate / genome_contig_file)

if (method=="sophisticated"):
    if(output_path_sophisticated.exists()==False):
        os.makedirs(output_path_sophisticated)
    score_threshold=sys.argv[4].split()[0]
    n_threshold=sys.argv[4].split()[1]
    fraction_threshold=sys.argv[4].split()[2]
    n_scc_threshold=sys.argv[4].split()[3]
    generate_output(genome_contig_file=genome_contig_file, score_threshold=score_threshold, n_threshold=n_threshold, fraction_threshold=fraction_threshold, n_scc_threshold=n_scc_threshold)
elif (method=="intermediate"):
    if(output_path_intermediate.exists()==False):
        os.makedirs(output_path_intermediate)
    score_threshold=sys.argv[4].split()[0]
    n_threshold=sys.argv[4].split()[1]
    initial_n_threshold=sys.argv[4].split()[2]
    distance_threshold=sys.argv[4].split()[3]
    generate_output_intermediate(genome_contig_file, score_threshold=score_threshold, n_threshold=n_threshold, initial_n_threshold=initial_n_threshold, distance_threshold=distance_threshold)
elif (method=="baseline"):
    if(output_path_baseline.exists()==False):
        os.makedirs(output_path_baseline)
    n_threshold=sys.argv[4]
    generate_output_baseline(genome_contig_file, n_threshold=n_threshold)
else: # simple
    if(output_path_simple.exists()==False):
        os.makedirs(output_path_simple)
    score_threshold=sys.argv[4].split()[0]
    #print(sys.argv[3])
    #print(sys.argv[4])
    n_threshold=sys.argv[4].split()[1]
    generate_output_simple(genome_contig_file, score_threshold=score_threshold, n_threshold=n_threshold)




