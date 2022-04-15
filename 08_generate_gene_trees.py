#!/usr/bin/env python
# coding: utf-8

"""Generates gene trees using neighbour joining (NJ), given MSAs (multiple sequence alignements) of all orthogroups, for all species.


For all species in the dataset, unrooted gene trees are generated using NJ, in newick format, for all orthogroups present in this species.

Input files: 
<orthogroup>.txt: this is a FASTA file containing the supermatrix, which is a concatenated alignment of the amino acid sequences of all core genes of all genomes belonging to a species.
filtered_pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output file:
<orthogroup>.nw: this is a treefile in newick format containing the gene tree of an orthogroup in a particular species.
"""

# import statements
import pandas as pd
from pathlib import Path
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# script parameters
project_path=Path().resolve().parent
msa_path=project_path / "results" / "intermediate" / "msas"
output_path=project_path / "results" / "intermediate" / "gene_trees"

# import data
pangenome=pd.read_csv(project_path / "results" / "intermediate" / "filtered_pangenome.csv")

def generate_genetree(orthogroup, species):
    """Generates the gene tree (using NJ, in newick format) belonging to a specific orthogroup in a specific species."""
    orthogroup_file=str(orthogroup+".txt")
    aln = AlignIO.read(open(msa_path / species / orthogroup_file), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    njtree = constructor.nj(dm)
    Phylo.write(njtree, output_path / species / str(orthogroup+".nw"), "newick")

# generate gene tree for all species, for all orthogroups within that species
for species in pangenome.gtdb_species.drop_duplicates():
    if(exists(output_path / species)==False):
        !mkdir $output_path/$species
        for orthogroup in pangenome.orthogroup.drop_duplicates():
            print(orthogroup)
            generate_genetree(orthogroup, species)
    print(str(species + " done"))





