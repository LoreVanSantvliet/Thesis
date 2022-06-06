#!/usr/bin/env python
# coding: utf-8

"""Calculates mobility scores of orthogroups.

Calculates mobility scores of orthogroups of a phylogenetic clade, based on the pangenome file of this clade and metadata regarding the genomes present in this clade. 

Input files: 
genomes_metadata.csv: this is a csv file containing the following columns: genome, gtdb_species, gtdb_genus, gtdb_family, quality, checkm_completeness, checkm_contamination, ncbi_isolation_source, ncbi_strain_identifiers, gtdb_representative,species. Only the first two columns are required.
pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output file:
mobility_frame.csv: this is a csv file containing the following columns: orthogroup, species_count, genomes_count, accessory_count.
"""

# import statements
import pandas as pd
from pathlib import Path
import numpy as np

# script parameters
genomes_per_species = 10 # minimal number of genomes per species
threshold_accessory = 0.9 # fraction of genomes per species an orthogroup can be present in at max. to be considered accessory
min_contigsize = 20 # minimal contig size (expressed in number of genes present in this contig) to take into account
project_path = Path().resolve().parent.parent
path_genomes = project_path / "data" / "genomes_metadata.csv"
path_pangenome = project_path / "data" / "pangenome.tsv"
output_path = project_path / "results" / "intermediate" / "mobility_frame.csv"

# import data
genomes = pd.read_csv(path_genomes).loc[:,'genome':'gtdb_species']
pangenome = pd.read_csv(path_pangenome, delimiter="\t", header=None)
pangenome.columns = ["gene", "genome", "orthogroup"]

def count_genomes_per_species(genomes):
    """Counts the number of genomes present per species present in a genomes metadata file."""
    counts = genomes.groupby(by="gtdb_species", as_index=False).count()
    counts.columns = ['gtdb_species', 'counts']
    return counts

def filter_species(genomes, filter_threshold):
    """Filters the species to contain a minimal number of genomes per species."""
    counts = count_genomes_per_species(genomes)
    genomes_filtered = pd.merge(genomes, counts[counts.counts >= filter_threshold], on='gtdb_species', how='right').sort_values('gtdb_species')
    return genomes_filtered

def count_orthogroups_in_species(df):
    """Counts the frequency of orthogroups in species (in how many species is this orthogroup present), which is an inverse mobility score."""
    orthocounts_species = df.loc[:,["orthogroup","gtdb_species"]].drop_duplicates().groupby(by="orthogroup", as_index=False).count()
    orthocounts_species.columns = ["orthogroup", "species_count"]
    return orthocounts_species

def count_orthogroups_in_genomes(df):
    """Counts the frequency of orthogroups in genomes (in how many genomes is this orthogroups present), which is an inverse mobility score."""
    orthocounts_genomes = df.loc[:,["orthogroup","genome"]].drop_duplicates().groupby(by="orthogroup", as_index=False).count()
    orthocounts_genomes.columns = ["orthogroup", "genomes_count"]
    return orthocounts_genomes

def count_orthogroups_in_genomes_per_species(df):
    """Counts the frequency of orthogroups per species (in how many genomes per species is this orthogroup present)."""
    orthocounts_genomes_per_species = df.loc[:,["orthogroup","gtdb_species","genome"]].drop_duplicates().groupby(by=["gtdb_species", "orthogroup"], as_index=False).count()
    orthocounts_genomes_per_species.columns = ["gtdb_species", "orthogroup", "orthocounts"]
    return orthocounts_genomes_per_species

def count_accessory_in_species(df):
    """Counts the frequency of orthogroups in the accessory genomes of species, which is a mobility score."""
    counts = count_genomes_per_species(genomes)
    orthocounts_genomes_per_species = count_orthogroups_in_genomes_per_species(df)
    orthocounts_accessory_per_species = pd.merge(orthocounts_genomes_per_species, counts, on='gtdb_species', how='left')
    orthocounts_accessory_per_species.columns = ["gtdb_species", "orthogroup", "orthocounts", "genome_counts"]
    orthocounts_accessory_per_species = orthocounts_accessory_per_species[orthocounts_accessory_per_species.orthocounts <= threshold_accessory * orthocounts_accessory_per_species.genome_counts]
    orthocounts_accessory_per_species = orthocounts_accessory_per_species.iloc[:, 0:2].groupby(by="orthogroup", as_index=False).count()
    orthocounts_accessory_per_species.columns = ["orthogroup", "accessory_count"] 
    orthocounts_accessory_per_species.accessory_count = orthocounts_accessory_per_species.accessory_count
    return orthocounts_accessory_per_species
   
def count_core_in_species(df):
    """Counts the frequency of orthogroups in the core genomes of species."""
    counts = count_genomes_per_species(genomes)
    orthocounts_genomes_per_species = count_orthogroups_in_genomes_per_species(df)
    orthocounts_core_per_species = pd.merge(orthocounts_genomes_per_species, counts, on='gtdb_species', how='left')
    orthocounts_core_per_species.columns = ["gtdb_species", "orthogroup", "orthocounts", "genome_counts"]
    orthocounts_core_per_species = orthocounts_core_per_species[orthocounts_core_per_species.orthocounts > threshold_accessory * orthocounts_core_per_species.genome_counts]
    orthocounts_core_per_species = orthocounts_core_per_species.iloc[:, 0:2].groupby(by="orthogroup", as_index=False).count()
    orthocounts_core_per_species.columns = ["orthogroup", "core_count"] 
    orthocounts_core_per_species.core_count = orthocounts_core_per_species.core_count
    return orthocounts_core_per_species
    
def generate_mobility_frame(df):
    """Generates a dataframe with the different mobility scores per orthogroup."""
    orthocounts_species = count_orthogroups_in_species(df)
    orthocounts_genomes = count_orthogroups_in_genomes(df)
    orthocounts_accessory = count_accessory_in_species(df)
    orthocounts_core = count_core_in_species(df)
    mobility_frame = pd.merge(orthocounts_species, orthocounts_genomes, on="orthogroup")
    mobility_frame = pd.merge(mobility_frame, orthocounts_accessory, on="orthogroup", how="left").fillna(0)
    mobility_frame = pd.merge(mobility_frame, orthocounts_core, on="orthogroup", how="left").fillna(0)
    mobility_frame['accessory_fraction']=mobility_frame.accessory_count.astype(float).divide(mobility_frame.species_count.astype(float))
    mobility_frame.iloc[:,1] = np.reciprocal(mobility_frame.iloc[:,1].astype(float))
    mobility_frame.iloc[:,2] = np.reciprocal(mobility_frame.iloc[:,2].astype(float))
    return(mobility_frame)

genomes_filtered = filter_species(genomes=genomes, filter_threshold=genomes_per_species)
full = pd.merge(genomes_filtered, pangenome, on='genome', how='left')
mobility_frame = generate_mobility_frame(full)
mobility_frame.to_csv(output_path)






