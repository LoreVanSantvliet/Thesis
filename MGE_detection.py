#!/usr/bin/env python
# coding: utf-8

"""Performs MGE detection in genomes.

Uses the phyletic distribution-based mobility scores of orthogroups to detect MGEs (mobile genetic elements) in contigs. The MGEs are returned in tabular format, where all gene numbers (as appearing in the contig of coice) belonging to a particular MGE are listed. The way this script sould be run is as follows: python MGE_detection.py <genome>-<contig>.csv, where the <genome>-<contig>.csv file is present in the input path (see script parameters).

Input files:
<genome>-<contig>.csv: this is a csv file containing the following columns: order (order of this particular gene in this contig), orthogroup (orthogroup ID), gene (gene ID), accessory_fraction (mobility score), count (how many times does this orthogroup appear in this genome), accessory (is this particular gene accessory in the species to which the genome belongs). 

Output files:
<genome>-<contig>.csv: this is a csv file containing the following columns: contig, MGE (integer indicating the MGE to which these particular gene belongs), gene_nr.
"""

# import statements
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# script parameters
score_threshold = 0.3 # threshold for a mobility score to be considered 'elevated'
n_threshold=10
fraction_threshold=0.1
project_path = Path().resolve().parent
input_path = project_path / "results" / "mobility_files"
output_path = project_path / "results" / "MGE_files_contig"

# import data
if(output_path.exists()==False):
    os.mkdir(output_path)
genome_contig_file=sys.argv[1]
frame=pd.read_csv(input_path / genome_contig_file, index_col=0)

def end_detection(mob_frame):
    # The index of the end boundary gene of the MGE equals the index of the last gene not belonging to the MGE (returned by applying the get_indices_lengths function on the mobility frame in reverse order) plus the length of the MGE. The gene_nr of thi bundary equals this index +1.
    genes,lengths=get_indices_lengths(mob_frame[::-1])
    end_genes=list(np.array(genes)+np.array(lengths))
    return [gene+1 for gene in end_genes][::-1]

def start_detection(mob_frame):
    # The index of the start boundary gene of the MGE equals the index of the first gene not belonging to the MGE (returned by applying the get_indices_lengths function on the mobility frame) minus the length of the MGE. The gene_nr of this bundary equals this index +1.
    genes,lengths=get_indices_lengths(mob_frame)
    start_genes=list(np.array(genes)-np.array(lengths))
    return [gene+1 for gene in start_genes]

def get_indices_lengths(frame, score_threshold=score_threshold, fraction_threshold=fraction_threshold, n_threshold=n_threshold):
    """Returns the index of the first gene not belonging to an MGE anymore when looping through consecutive genes in a contig, and length of the MGE."""
    i=0
    exc=0
    indices=[]
    lengths=[]
    for index,row in frame.iterrows():
        consecutive=(row.accessory_fraction>score_threshold)
        if consecutive:
            # look for consecutive genes with an elevated mobility score
            i=i+1
        else:
            exc=exc+1
            #exc are exceptions, so genes that do not have an elevated mobility score
            if(row.accessory|(row['count']>1)|(exc<fraction_threshold*i)):
                # in case these exceptions are multi-copy or accessory, or when they do not appear too close to the start of the potential MGE, they are tolerated
                i=i+1
            else:
                # otherwise, they are not
                if i>n_threshold:
                    # if the number of consecutive genes that belong to a potential MGE is large enough, the potential MGE is predicted to be a true MGE
                    # in this case, the index and length of the MGE are saved
                    indices.append(index)
                    lengths.append(i)
                # parameters are set to zero again, to find a new (potential) MGE
                i=0
                exc=0
    return indices,lengths

def generate_output(genome_contig_file):
    """Generates csv files, containing a columns with the contig, a column with a number indicating an individual MGE, and a columns with gene_nrs, indicating the genes that belong to a particular MGE."""
    contig_extension=genome_contig_file.split('-')[1]
    contig=str(contig_extension.split('.')[0]+'.'+contig_extension.split('.')[1])
    MGE_list=[]
    mob_frame=pd.read_csv(input_path / genome_contig_file, index_col=0)
    start_genes=start_detection(mob_frame)
    end_genes=end_detection(mob_frame)
    if (len(start_genes)>0)&(len(end_genes)>0):
        if(len(start_genes)!=len(end_genes)):
            print("PROBLEM: start and end boundary mismatch")
        for i in range(0,len(start_genes)):
            for gene in range(start_genes[i], end_genes[i]+1):
                MGE_list.append({'contig':contig,'MGE':i+1,'gene_nr':gene})
    MGE_frame = pd.DataFrame.from_records(MGE_list)
    if(not MGE_frame.empty):
        MGE_frame.to_csv(output_path / genome_contig_file)       

generate_output(genome_contig_file)

