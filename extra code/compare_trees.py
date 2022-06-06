#!/usr/bin/env python
# coding: utf-8

"""Compares gene and species trees using treeKO distances."""

# import statements
#import pandas as pd
from pathlib import Path
from ete3 import PhyloTree
#from ete3 import Tree
import os
import sys
from csv import writer
import subprocess
import numpy

numpy.seterr(all='raise')

# script parameters
project_path = Path().resolve().parent.parent
path_gene_trees = project_path / "results" / "intermediate" / "gene_trees_mod"
path_gene_trees_alt = project_path / "results" / "intermediate" / "gene_trees_mod_alt"
path_species_trees = project_path / "results" / "intermediate" / "species_trees_mod"
output_path = project_path / "results" / "intermediate" / "tree_mobility_frames"

species_1=sys.argv[1]
species_2=sys.argv[2]
orthogroup_file=sys.argv[3]
orthogroup=orthogroup_file.split('.')[0]
species=str(species_1+" "+species_2)
#orthogroup_file=str(orthogroup+".nw")
species_input_file=str(species+".nw")
species_output_file=str(species+".csv")

def remove_duplications_in_endnodes(gene_tree):
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
#print("Newick: " + str(path_gene_trees / species / orthogroup_file))
#print("Newick: " + str("'" + str(path_gene_trees / species / orthogroup_file) + "'"))
gene_tree=PhyloTree(gene_tree)
gene_tree.set_species_naming_function(lambda node: node.name.split('-')[0])
species_tree=PhyloTree(species_tree)
species_tree.set_species_naming_function(lambda node: node.name)
print(orthogroup_file)



#print(gene_tree)
#print(species_tree)
#print(gene_tree.get_species())
#print(species_tree.get_species())
#print(gene_tree.write())
#print(species_tree.write())
#if(len(gene_tree.split_by_dups())==1):
#    path=str("../results/intermediate/gene_trees_mod/" + species_1 + "\ " + species_2 + "/" + orthogroup_file)
#    path=path_gene_trees / species / orthogroup_file
#    subprocess.call(["sed", "-i", "s/-[A-Z]*[0-9]*\.[0-9]*_[0-9]*//g", path])
#    with open(path_gene_trees / species / orthogroup_file, 'r') as f:
#        gene_tree=f.read()
#        gene_tree=PhyloTree(gene_tree)
#        gene_tree.set_species_naming_function(lambda node: node.name)
    #row=[species, orthogroup]
    # subprocess.call(["./modify_nodes.py", argument1], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

    #with open(project_path / "results" / "intermediate" / tree_mod_file, 'a') as f_object:
    #    writer_object = writer(f_object)
    #    writer_object.writerow(row)
    #    f_object.close()

# hypothesis: trees of one node cause mean of empty slice error

#with open('first_try.csv', 'w') as doc:
#    doc.write("species,orthogroup \n")
#with open('second_try.csv', 'w') as doc:
#    doc.write("species,orthogroup \n")
#with open('third_try.csv', 'w') as doc:
#    doc.write("species,orthogroup \n")
#with open('failed.csv', 'w') as doc:
#    doc.write("species,orthogroup \n")

#f_object = open('first_try.csv', 'a')


try:
    comp=gene_tree.compare(species_tree, has_duplications=True, unrooted=True)
    #print(gene_tree)
    #print(species_tree) 
    mob_score=comp.get('treeko_dist')
    #print("reached filewriting 1")
    with open('first_try.csv', 'a') as f_object:
        writer_object = writer(f_object)
        writer_object.writerow([species, orthogroup])
    #f_object.write(str(species + ", " +  orthogroup + "\n"))
    #print("closing")
except BaseException:
    try:
        remove_duplications_in_endnodes(gene_tree)
        comp=gene_tree.compare(species_tree, has_duplications=True, unrooted=True)
        mob_score=comp.get('treeko_dist')
        print("reached filewriting 2")
        #f_object = open('second_try.csv', 'a')
            #writer_object = writer(f_object)
                #writer_object.writerow([species, orthogroup])
        #f_object.write(str(species + ", " +  orthogroup + "\n"))
        #f_object.close()
        with open('second_try.csv', 'a') as f_object:
            writer_object = writer(f_object)
            writer_object.writerow([species, orthogroup])
            #f_object.close()
    except BaseException:
        try:
            comp=gene_tree.compare(species_tree, has_duplications=False, unrooted=True)
            mob_score=comp.get('norm_rf')
            print("reached filewriting 3")
            #f_object = open('first_try.csv', 'a')
                #writer_object = writer(f_object)
                    #writer_object.writerow([species, orthogroup])
            #f_object.write(str(species + ", " +  orthogroup + "\n"))
            #f_object.close()

            with open('third_try.csv', 'a') as f_object:
                writer_object = writer(f_object)
                writer_object.writerow([species, orthogroup])
                #f_object.close()
        except BaseException:
            print("reached filewriting 4")
            #f_object = open('first_try.csv', 'a')
            #    writer_object = writer(f_object)
            #    writer_object.writerow([species, orthogroup])
            #f_object.write(str(species + ", " +  orthogroup + "\n"))
            #f_object.close()

            with open('failed.csv', 'a') as f_object:
                writer_object = writer(f_object)
                writer_object.writerow([species, orthogroup])
                #f_object.close()

    #print(gene_tree)
    #print(species_tree)

    #with open(path_gene_trees_alt / species / orthogroup_file, 'r') as f:
    #    gene_tree_alt=f.read()
    #gene_tree_alt=PhyloTree(gene_tree_alt)
    #gene_tree_alt.set_species_naming_function(lambda node: node.name.split('-')[0])
    #print(species_tree)
    #print(gene_tree_alt)
    #print(gene_tree_alt.write())
    #comp=gene_tree_alt.compare(species_tree, has_duplications=True, unrooted=True)
print(mob_score)
row=[orthogroup, mob_score]
with open(output_path/species_output_file, 'a') as f_object:
    writer_object = writer(f_object)
    writer_object.writerow(row)
    #f_object.close()
print(orthogroup_file + " finished")


f_object.close()
