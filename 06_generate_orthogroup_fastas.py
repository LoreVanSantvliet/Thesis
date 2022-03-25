#!/usr/bin/env python
# coding: utf-8

"""Generates FASTA files of orthogroups, per species.

Collects nulceotide sequences from ffn files, where they are ordered per genome, and transfers them to FASTA files, where they are ordered per orthogroup. These FASTA files are generated per orthigroup, per species.

Input files: 
<species>.fasta: this is a FASTA file containing the supermatrix, which is a concatenated alignment of the amino acid sequences of all core genes of all genomes belonging to a species.
filtered_pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output file:
<species>.nw: this is a treefile in newick format containing the species tree of a species.
"""

# import statements
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# script parameters
project_path=Path().resolve().parent
genomes_path=project_path / "results" / "intermediate" / "filtered_ffns"
orthogroup_seq_path=project_path / "results" / "intermediate" / "gene_tree_sequences"

# import data
pangenome=pd.read_csv(project_path / "results" / "intermediate" / "filtered_pangenome.csv")


def read_ffn(genome):
    """Reads ffn files and generates a dataframe with columns 'gene' and 'sequence', containing all genes and sequences that are present in the ffn file of a genome of choice."""
    genome_file=genome+".ffn"
    df = pd.read_table(
    genomes_path / genome_file,
    engine = 'c',
    lineterminator = '>',
    skiprows =1,
    names = ['raw']
    )

    # The first line break ('\n') separates Column 0 from Column 1
    df[['header','sequence']] = pd.DataFrame.from_records(df.raw.apply(lambda s: s.split(maxsplit=1, sep="\n")))

    # All subsequent line breaks (which got left in Column 1) should be ignored
    df['sequence'] = df['sequence'].apply(lambda s: s.replace('\n',''))
    
    # Only the first word of the header is relevant.
    df['gene'] = pd.DataFrame.from_records(df.header.apply(lambda s: s.split(maxsplit=1)))[0]
    return df[['gene', 'sequence']]

def add_to_file(orthogroup, gene, genome, sequence, species):
    """Adds a sequence of a gene in a particular genome, belonging to a species to the correct orthogroup FASTA file. The header of this sequence in the newly generated FASTA file is the genome, followed by a dash and the gene name."""
    file_title = orthogroup_seq_path / species / str(orthogroup+".txt")
    file_object = open(file_title, 'a')
    file_object.write(str(">"+genome+"-"+gene+"\n"))
    file_object.write(str(sequence+"\n"))
    file_object.close()


def add_files(genome, species):
    """Appends all nucleotide sequences belonging to a genome, to the different orthogroup FASTA files of the correct species."""
    df=read_ffn(genome)
    length=len(df)
    for j in range(0, df.shape[0]):
        add_to_file(pangenome.orthogroup.iloc[j], pangenome.gene.iloc[j], genome, df.sequence.iloc[j], species)

# For all genomes, for all species, generate the orthogroup FASTA files.
for species in pangenome.gtdb_species.drop_duplicates():
    for genome in pangenome.genome[pangenome.gtdb_species == species].drop_duplicates():
        add_files(genome, species)
        print(str(genome + " done!"))
    print(str(species + " done"))
    