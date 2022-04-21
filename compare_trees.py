#!/usr/bin/env python
# coding: utf-8

"""Visualizes mobility scores of orthogroups over a genome.

Visualizes the mobility score of orthogroups  over a genome of choice to detect mobile genetic elements (MGEs).

Input files: 
genomes_metadata.csv: this is a csv file containing the following columns: genome, gtdb_species, gtdb_genus, gtdb_family, quality, checkm_completeness, checkm_contamination, ncbi_isolation_source, ncbi_strain_identifiers, gtdb_representative,species. Only the first two columns are required.
mobility_frame.csv: this is a csv file containing the following columns: orthogroup, species_count, genomes_count, accessory_fraction.
pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output files: 
<genome>-<contig>.png: mobility plot showing the mobility scores for cosecutive genes of a genome of choice.
<genome>-<contig>.csv: mobility frame for consecutive genes of a genome of choice.
"""

# import statements
import pandas as pd
from pathlib import Path
from ete3 import PhyloTree
import os

# script parameters
project_path = Path().resolve().parent
path_gene_trees = project_path / "results" / "intermediate" / "gene_trees_mod"
path_species_trees = project_path / "results" / "intermediate" / "species_trees"
output_path(project_path / "results" / "tree_mobility_frames")

orthogroup=sys.argv[1]
species_1=sys.argv[2]
species_2=sys.argv[3]
species=str(species_1+" "+species_2)
orthogroup_file=str(orthogroup+".nw")
species_input_file=str(species+".nw")
species_output_file=str(species+".csv")

gene_tree=PhyloTree(path_gene_trees / orthogroup_file)
gene_tree.set_species_naming_function(lambda node: node.name)
species_tree=PhyloTree(path_species_trees / species_input_file)
species_tree.set_species_naming_function(lambda node: node.name)
comp=gene_tree.compare(species_tree, has_duplications=True)
mob_score=comp.get('treeko_dist')
row=[species, mob_score]
with open(output_path/species_output_file, 'a') as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(row)
    f_object.close()






