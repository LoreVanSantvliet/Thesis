#!/usr/bin/env python
# coding: utf-8

"""Generates species trees using neighbour joining (NJ), given supermatrices, which are concatenated alignments of core genes.

For all species in the dataset, a distancematrix is created (biopython), from which an unrooted species tree is generated using NJ (biopython), in newick format.

Input files: 
<species>.fasta: this is a FASTA file containing the supermatrix, which is a concatenated alignment of the amino acid sequences of all core genes of all genomes belonging to a species.
filtered_pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output file:
<species>.nw: this is a treefile in newick format containing the species tree of a species.
"""

# import statements
import pandas as pd
from pathlib import Path
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# script parameters
project_path=Path().resolve().parent.parent
msa_path=project_path / "results" / "intermediate" / "supermatrix"
output_path=project_path / "results" / "intermediate" / "species_trees"

# import data
pangenome=pd.read_csv(project_path / "results" / "intermediate" / "filtered_pangenome.csv")

def generate_speciestree(species):
    """Generates a species tree, given the name of a species, by navigating to the correct location, to read in the supermatrix of this species, and uses NJ (neighbour joining) to generate an unrooted species tree in newick format."""
    species_file=str(species+".fasta")
    aln = AlignIO.read(open(msa_path / species_file), 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    njtree = constructor.nj(dm)
    Phylo.write(njtree, output_path / str(species+".nw"), "newick")

# loop over all species and generate species trees
for species in pangenome.gtdb_species.drop_duplicates():
    print(species)
    generate_speciestree(species)





