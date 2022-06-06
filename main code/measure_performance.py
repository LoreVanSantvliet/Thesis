#!/usr/bin/env python
# coding: utf-8

"""Measures performance measures for a set of mges detected in a genome, when compared to a set of MGEs found by a tool.

Measures the performance (coverage, dispersion, recall, precision and F-score) of a set of mges detected, when compared to a set of mges found by the tools CONJscan and Phaster, as well as the joined results of these tools. Writes the results to specified files. To be run as follows: "python measure_performance.py self_path conjscan_path conjscan_file_path phaster_path phaster_file_path joined_file_path genome n_threshold"

Input files: 
<species>.csv (in tree_mobility frames): this is a csv file with columns orthogroup and tree_score. The tree_score column represents the mobility score of an orthogroup within one particular species.
genomes_species.csv: this is a csv file with columns gtdb_species and genome, indicating which genomes belong to which species. 
mobility_frame.csv: this is a csv file with the following columns: orthogroup, species_count, genomes_count, accessory_count, core_count, accessory_fraction. Accessory_fraction represents the phyletic distribution pattern-based mobility score.

Output files:
conjscan_performance.csv: this is a csv file with columns genome, n_threshold, coverage, dispersion, recall, precision and F-score, where a row represents a set of performance measures for MGE detection in a specified genome, with specified threshold values.
phaster_performance.csv: this is a csv file with columns genome, n_threshold, coverage, dispersion, recall, precision and F-score, where a row represents a set of performance measures for MGE detection in a specified genome, with specified threshold values.
joined_performance.csv: this is a csv file with columns genome, n_threshold, coverage, dispersion, recall, precision and F-score, where a row represents a set of performance measures for MGE detection in a specified genome, with specified threshold values.
"""

# Measures the performance (coverage, dispersion, recall, precision and F-score) of a set of mges detected, when compared to a set of mges found by a tool. Writes the results to a specified file. To be run as follows: "python measure_performance.py self_path conjscan_path conjscan_file_path phaster_path phaster_file_path joined_file_path genome n_threshold"

# import statements
import pandas as pd
from pathlib import Path
from csv import writer
import sys
import numpy as np
import os
import re


print(sys.argv)

# read input
method = sys.argv[1]
self_path = sys.argv[2]
conjscan_path = sys.argv[3]
conjscan_file_path = sys.argv[4]
phaster_path = sys.argv[5]
phaster_file_path = sys.argv[6]
joined_file_path = sys.argv[7]
genome = sys.argv[8]

print("Self_path: " + self_path)

if method == "baseline":
    n_threshold = sys.argv[9]
elif method == "simple":
    score_threshold = sys.argv[9].split()[0]
    n_threshold = sys.argv[9].split()[1]
elif method == "intermediate":
    score_threshold = sys.argv[9].split()[0]
    n_threshold = sys.argv[9].split()[1]
    initial_n_threshold = sys.argv[9].split()[2]
    distance_threshold = sys.argv[9].split()[3]
else:
    score_threshold = sys.argv[9].split()[0]
    n_threshold = sys.argv[9].split()[1]
    fraction_threshold = sys.argv[9].split()[2]
    n_scc_threshold = sys.argv[9].split()[3]

def get_performance_measures(tool_mges, self_mges):
    tool_MGEs=tool_mges.MGE.drop_duplicates()
    n_tool = len(tool_mges.loc[:,['contig', 'MGE']].drop_duplicates())
    tool_mges_detected=[]
    coverages_denominator=0
    coverages_numerator=0
    dispersions=[]
    for tool_MGE in tool_MGEs:
        detected_flag = False
        coverage_nr = 0
        dispersion_set=set()
        for tool_row in tool_mges.loc[tool_mges.MGE == tool_MGE].iterrows():
            #print("Error: " + tool_mges.contig.loc[tool_mges.MGE == tool_MGE])
            gene = tool_row[1][3]
            tool_contig = tool_row[1][1]
            for self_row in self_mges.iterrows():
                if (self_row[1][3] == gene) & (self_row[1][1] == tool_contig):
                    detected_flag = True
                    coverage_nr = coverage_nr + 1
                    dispersion_set.add(self_row[1][2])
        if(detected_flag):
            tool_mges_detected.append(tool_MGE)
            coverages_denominator=coverages_denominator+len(tool_mges.gene_nr.loc[tool_mges.MGE == tool_MGE])
            coverages_numerator=coverages_numerator+coverage_nr
            dispersions.append(dispersion_set)
    number_mges = []
    for mge_set in dispersions:
        number_mges.append(len(mge_set))
    average_dispersion=np.mean(number_mges)
    #print(tool_mges_detected)
    #print(coverages)
    #print(dispersions)     
    
    if not tool_mges_detected:
        print("No MGEs detected")
        if method=="intermediate":
            return [genome,score_threshold,n_threshold,initial_n_threshold,distance_threshold,n_tool,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        elif method=="baseline":
            return [genome,n_threshold,n_tool,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        elif method=="simple":
            return [genome,score_threshold,n_threshold,n_tool,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        else:
            return [genome,score_threshold, n_threshold, fraction_threshold, n_scc_threshold,n_tool,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
    
    else:
        if method=="intermediate":
            return [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold, n_tool, len(tool_mges_detected), len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        elif method=="baseline":
            return [genome, n_threshold, n_tool, len(tool_mges_detected), len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        elif method=="simple":
             return [genome, score_threshold, n_threshold, n_tool, len(tool_mges_detected), len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]
        else:
            return [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold, n_tool, len(tool_mges_detected), len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),average_dispersion, coverages_denominator, coverages_numerator]

def write_performance_line(self_path, tool_path, file_path):
    if (os.path.exists(tool_path)):
        tool_mges = pd.read_csv(tool_path)
        if (os.path.exists(self_path)):
            self_mges = pd.read_csv(self_path)
            row = get_performance_measures(tool_mges, self_mges)
        else:
            if method=="intermediate":
                row = [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold,len(tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
            elif method=="baseline":
                row = [genome, n_threshold,len(tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
            elif method=="simple":
                row = [genome, score_threshold, n_threshold,len(tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
            else:
                row = [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold,len(tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]

    else:
        print(tool_path + " does not exist")
        if (os.path.exists(self_path)):
            self_mges = pd.read_csv(self_path)
            if method=="intermediate":
                row = [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            elif method=="baseline":
                row = [genome, n_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            elif method=="simple":
                row = [genome, score_threshold, n_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            else:
                row = [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            i=0
        else:
            print(self_path + " does not exist")
            if method=="intermediate":
                row = [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold,0,0,0,None,0,0]
            elif method=="baseline":
                row = [genome, n_threshold,0,0,0,None,0,0]
            elif method=="simple":
                row = [genome, score_threshold, n_threshold,0,0,0,None,0,0]
            else:
                row = [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold,0,0,0,None,0,0]
    with open(file_path, 'a') as f_object:
        i=0
        for element in row:
            f_object.write(str(element))
            i=i+1
            if (i<len(row)):
                f_object.write(",")
        f_object.write("\n")


def write_joined_performance_line(self_path, conjscan_path, phaster_path, joined_file_path):
    if(os.path.exists(phaster_path)):
        phaster_mges = pd.read_csv(phaster_path, index_col=0)
    if(os.path.exists(conjscan_path)):
        conjscan_mges = pd.read_csv(conjscan_path, index_col=0)
    if(os.path.exists(phaster_path) & os.path.exists(conjscan_path)):
        
        #phaster_mges['gene_nr']=phaster_mges['gene_nr'].astype('float64')
        #print(phaster_mges['gene_nr'].dtype)	
        #print(phaster_mges['contig'].dtype)
        #conjscan_mges['gene_nr']=conjscan_mges['gene_nr'].astype('float64')
        #conjscan_mges['contig']=conjscan_mges['gene_nr'].astype(str)
        #print(conjscan_mges['gene_nr'].dtype)	
        #print(conjscan_mges['contig'].dtype)

        joined_tool_mges = pd.merge(phaster_mges, conjscan_mges, on = ['contig', 'gene_nr'], how='outer')
        joined_tool_mges['MGE']=0
        for index,row in joined_tool_mges.iterrows():
            if (np.isnan(row.MGE_x) == False):
                joined_tool_mges.loc[index,'MGE'] = row.MGE_x
            else:
                joined_tool_mges.loc[index,'MGE'] = row.MGE_y + np.max(joined_tool_mges.MGE_x)
    elif (os.path.exists(phaster_path)):
        joined_tool_mges = phaster_mges
    elif (os.path.exists(phaster_path)):
        joined_tool_mges = conjscan_mges
    else:
        joined_tool_mges=pd.DataFrame(columns=['index', 'contig', 'MGE', 'gene_nr'])    
    
    if (os.path.exists(self_path)):
        joined_tool_mges['index']=joined_tool_mges.index
        self_mges = pd.read_csv(self_path)
        row = get_performance_measures(joined_tool_mges.loc[:,['index', 'contig', 'MGE', 'gene_nr']], self_mges)
    
        if((os.path.exists(phaster_path)==False) & (os.path.exists(conjscan_path)==False)):
            if method=="intermediate":
                row = [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            elif method=="baseline":
                row = [genome, n_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            elif method=="simple":
                row = [genome, score_threshold, n_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
            else:
                row = [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold,0,0,len(self_mges.loc[:,['contig', 'MGE']].drop_duplicates()),None,0,0]
	
    else:
        if method=="intermediate":
            row = [genome, score_threshold,n_threshold,initial_n_threshold,distance_threshold,len(joined_tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
        elif method=="baseline":
            row = [genome, n_threshold,len(joined_tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
        elif method=="simple":
            row = [genome, score_threshold, n_threshold,len(joined_tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
        else:
            row = [genome, score_threshold, n_threshold, fraction_threshold, n_scc_threshold,len(joined_tool_mges.loc[:,['contig', 'MGE']].drop_duplicates()),0,0,None,0,0]
    
    #self_mges['index'] = self_mges.index
    #print(joined_tool_mges.loc[:,['index', 'contig', 'MGE', 'gene_nr']])
   
    with open(joined_file_path, 'a') as f_object:
        i=0
        for element in row:
            f_object.write(str(element))
            i=i+1
            if (i<len(row)):
                f_object.write(",")
        f_object.write("\n")
        
print("Conjscan analysis")
write_performance_line(self_path, conjscan_path, conjscan_file_path)
print("Phaster analysis")
write_performance_line(self_path, phaster_path, phaster_file_path)
print("Joined performance analysis")
write_joined_performance_line(self_path, conjscan_path, phaster_path, joined_file_path)



            
