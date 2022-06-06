#!/usr/bin/env python
# coding: utf-8


#import statements
import pandas as pd
from pathlib import Path
import sys

directory=sys.argv[1]

# script parameters
project_path=Path().resolve().parent.parent.parent
pangenome_path=project_path / "results" / "intermediate" / "filtered_pangenome.csv"
genomes_path=project_path / "results" / "intermediate" / directory / "genomes.csv"
result_path=project_path / "results" / "intermediate" / directory / "conjscan"

genomes = pd.read_csv(genomes_path, header=None)
pangenome = pd.read_csv(pangenome_path)
pangenome['order'] = pd.to_numeric(pangenome["gene"].str.split("_", n=1, expand=True)[1])
pangenome['contig'] = pangenome["gene"].str.split("_", n=1, expand=True)[0]

for genome_file in genomes.iloc[:,0]:
    print(genome_file)
    genome = str(genome_file.split('.')[0] + "." + genome_file.split('.')[1])
    pangenome_sel = pangenome[pangenome.genome==genome].sort_values(['contig', 'order']).reset_index()
    pangenome_sel.to_csv(result_path / str("order_" + genome_file))



