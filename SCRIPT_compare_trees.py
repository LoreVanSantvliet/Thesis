#!/usr/bin/env python
# coding: utf-8

"""Compares a gene tree with a species tree.

Compares a gene tree with a species tree, and returnsa distace measure which represents the difference between both trees, using ETE.

Input files: 
<species>.nw: this is a newick file containing the species tree.
<orthogroup>.nw (in species directory): this is a newick file containing the gene tree for a specific orthogroup on the species level.

Output files: 
<species>.csv: this is a csv file containing the mobility scores (treeko distances) for orthogroups, calculated in a particular species.
"""

# import statements
from pathlib import Path
from ete3 import PhyloTree
import os
import sys
from csv import writer
import numpy

numpy.seterr(all='raise')

# script parameters
project_path = Path().resolve().parent
path_gene_trees = project_path / "results" / "intermediate" / "gene_trees_mod"
path_gene_trees_alt = project_path / "results" / "intermediate" / "gene_trees_mod_alt"
path_species_trees = project_path / "results" / "intermediate" / "species_trees_mod"
output_path = project_path / "results" / "intermediate" / "tree_mobility_frames"

# read data
species_1=sys.argv[1]
species_2=sys.argv[2]
orthogroup_file=sys.argv[3]
orthogroup=orthogroup_file.split('.')[0]
species=str(species_1+" "+species_2)
species_input_file=str(species+".nw")
species_output_file=str(species+".csv")

def remove_duplications_in_endnodes(gene_tree):
    """Removes one of both sister endnodes in case they are from the same genome, and thus are the result of a gene duplication."""
    seen = set()
    duplicates = [name for name in gene_tree.get_leaf_names() if name in seen or seen.add(name)]
    for name in duplicates:
        node_collection=gene_tree.search_nodes(name=name)
        num=len(node_collection)
        check_until=num-1
        for i in range(0, num):
            for j in range(i+1, num):
                sis=node_collection[i].get_sisters()
                if(sis):
                    if node_collection[i].get_sisters()[0]==node_collection[j]:
                        node_collection[i].delete()

print(orthogroup_file + " started")

with open(path_gene_trees / species / orthogroup_file, 'r') as f:
    gene_tree=f.read()
with open(path_species_trees / species_input_file, 'r') as f:
    species_tree=f.read()
gene_tree=PhyloTree(gene_tree)
gene_tree.set_species_naming_function(lambda node: node.name.split('-')[0])
species_tree=PhyloTree(species_tree)
species_tree.set_species_naming_function(lambda node: node.name)
print(orthogroup_file)

try:
    comp=gene_tree.compare(species_tree, has_duplications=True, unrooted=True)
    mob_score=comp.get('treeko_dist')
except BaseException:
    try:
        remove_duplications_in_endnodes(gene_tree)
        comp=gene_tree.compare(species_tree, has_duplications=True, unrooted=True)
        mob_score=comp.get('treeko_dist')
    except BaseException:
        comp=gene_tree.compare(species_tree, has_duplications=False, unrooted=True)
        mob_score=comp.get('norm_rf')
print(mob_score)
row=[orthogroup, mob_score]
with open(output_path/species_output_file, 'a') as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(row)
    f_object.close()
print(orthogroup_file + " finished")

