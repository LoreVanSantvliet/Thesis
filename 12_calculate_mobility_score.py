#!/usr/bin/env python
# coding: utf-8

"""Calculates tree based mobility scores of orthogroups.

Calculates mobility scores of orthogroups of a phylogenetic clade, based on the weighted averages of treeko distances between species and gene trees for these orthogroups in the different species. For these weighted averages, the weight given to a treeko distance in a particular species is proportional to the number of genomes present in this species.

Input files: 
<species>.csv (in tree_mobility frames): this is a csv file with columns orthogroup and tree_score. The tree_score column represents the mobility score of an orthogroup within one particular species.
genomes_species.csv: this is a csv file with columns gtdb_species and genome, indicating which genomes belong to which species. 
mobility_frame.csv: this is a csv file with the following columns: orthogroup, species_count, genomes_count, accessory_count, core_count, accessory_fraction. Accessory_fraction represents the phyletic distribution pattern-based mobility score.

Output file:
mobility_frame.csv:this is a csv file like the input file names 'mobility_frame.cv', with extra columns counts_mob and tree_score, indicating the number of genomes belonging to a species where this orthogroup is present, and the tree-based mobility score, respectively.
"""

# import statements
import pandas as pd
from pathlib import Path
import numpy as np

# script parameters
project_path = Path().resolve().parent
input_path=project_path / "results" / "intermediate"
path_mobility = input_path / "tree_mobility_frames"
path_genomes = input_path / "genomes_species.csv"
output_path = project_path / "results" / "intermediate" / "mobility_frame.csv"

# import data
mob_frame=pd.read_csv(output_path)
genomes_frame=pd.read_csv(path_genomes)
counts = genomes_frame.groupby(by = "gtdb_species", as_index = False).count()
counts.columns = ['gtdb_species', 'counts']
counts['mob']=0
mob_frame['counts_mob']=0
mob_frame['tree_score']=0

def calculate_mob_score(counts, path_mobility, mob_frame):
    """Calculates the tree-based mobility score of orthogroups by taking the (genome-)weighted average of treeko distances for all species.
    
    args:
    counts: dataframe containing the columns 'gtdb_species' and 'counts', which indicates the number of genomes for the species.
    path_mobility: path specifying the location of the mobility csv files for all species.
    mob_frame: frame containing the phyletic dicstribution-based mobility scores for all orthogroups. To this frame, the tree-based mobility scores will be added in a new column.
    """
    for index,row in counts.iterrows():
        species=row.gtdb_species
        print(species)
        mob_file=pd.read_csv(path_mobility / str(species+".csv"))
        for orthogroup in mob_file.orthogroup:
            if (np.isnan(mob_frame.loc[mob_frame.orthogroup==orthogroup, 'tree_score'].item())==False):
                mob_frame.loc[mob_frame.orthogroup==orthogroup, 'counts_mob']=mob_frame.counts_mob[mob_frame.orthogroup==orthogroup] + counts.counts[index]
                mob_frame.loc[mob_frame.orthogroup==orthogroup, 'tree_score']=mob_frame.loc[mob_frame.orthogroup==orthogroup, 'tree_score'].item() + (counts.counts[index])*(mob_file.treeko_dist[mob_file.orthogroup==orthogroup].item())
        print(species + " done")
    mob_frame.tree_score=mob_frame.tree_score/mob_frame.counts_mob
    return mob_frame

mob_frame=calculate_mob_score(counts, path_mobility, mob_frame)
mob_frame.to_csv(output_path)
