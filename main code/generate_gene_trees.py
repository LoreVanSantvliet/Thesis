#!/usr/bin/env python
# coding: utf-8

"""Generates an unrooted gene tree using neighbour joining (NJ), given an MSA (multiple sequence alignment) for a user supplied orthogroup and species.

A distance matrice is calculated (biopython) from the MSA. An unrooted gene tree in newick format is generated from this distance matrix using NJ (biopython).

Input files: 
<orthogroup>.txt: this is a FASTA file containing the supermatrix, which is a concatenated alignment of the amino acid sequences of all core genes of all genomes belonging to a species.

Output file:
<orthogroup>.nw: this is a treefile in newick format containing the gene tree of an orthogroup in a particular species.
"""

# import statements
from pathlib import Path
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import os
import sys

# script parameters
project_path=Path().resolve().parent.parent
msa_path=project_path / "results" / "intermediate" / "msas"
output_path=project_path / "results" / "intermediate" / "gene_trees"

def generate_genetree(orthogroup, species, species_2):
    """Generates the gene tree (using NJ, in newick format) belonging to a specific orthogroup in a specific species."""
    species_string=str(species + " " + species_2)
    orthogroup_string=orthogroup.split(sep=".")[0]

    if(os.path.exists(output_path / species_string)==False):
        os.mkdir(output_path / species_string)

    if(os.path.exists(output_path / species_string / orthogroup)==False):
        aln = AlignIO.read(open(msa_path / species_string / orthogroup), 'fasta')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        njtree = constructor.nj(dm)
        Phylo.write(njtree, output_path / species_string / str(orthogroup_string+".nw"), "newick")

generate_genetree(sys.argv[1], sys.argv[2], sys.argv[3])



