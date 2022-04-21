#!/usr/bin/env python
# coding: utf-8

"""Parses a ConjScan output file.

Parses a user supplied, raw ConjScan report into a file containing contigs and gene numbers of detected conjugation systems. This script should be run as follows: "python ParseSconjScan.py <genome>.txt".

Input files: 
<genome>.txt: this is a raw ConjScan report.
Can be obtained from https://galaxy.pasteur.fr/root?tool_id=toolshed.pasteur.fr%2Frepos%2Fodoppelt%2Fconjscan%2FConjScan%2F1.0.2, where genomic .faa files (translated CDS) are uploaded, and the following options need to be chosen:
The type of dataset to deal with: ordered replicon
Conjugative element to detect: typeFATA
Tune or leave default values to Hmmer options: default

Output file:
<genome>.txt: this is a parsed ConjScan report, in csv format, containing columns contig and gene.
"""

# import statements
import pandas as pd
from pathlib import Path
import sys

# script parameters
project_path=Path().resolve().parent.parent.parent
input_path=project_path / "results" / "intermediate" / "benchmarking" / "conjscan_raw"
result_path=project_path / "results" / "intermediate" / "benchmarking" / "conjscan_parsed"

genome_file=sys.argv[1]

file = pd.read_csv(input_path / genome_file, sep="\t")
file["contig"]=pd.DataFrame.from_records(file['hit_id'].apply(lambda s: s.split(sep='|')))[1]
file["contig"]=pd.DataFrame.from_records(file['contig'].apply(lambda s: s.split(sep='_')))[0] + "_" +pd.DataFrame.from_records(file['contig'].apply(lambda s: s.split(sep='_')))[1]
file["MGE"]=pd.DataFrame.from_records(file['sys_id'].apply(lambda s: s.split(sep='_')))[3]
file["gene"]=file.hit_pos
file_final = file.loc[:,['contig', 'MGE', 'gene']]
file_final.to_csv(result_path / str(genome_file))
