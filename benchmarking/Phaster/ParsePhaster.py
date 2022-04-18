#!/usr/bin/env python
# coding: utf-8

"""Parses a Phaster output file.

Parses a user-supplied, Phaster output file into a file containing contigs and gene numbers of detceted prophages. This script should be run as follows: "python ParsePhaster.py <genome>.txt"

Input files: 
<genome>.txt: this is a raw Phasetr output file, where newlines are fixed and the '----' line is replaced (see fix_newline.sh script).

Output file:
<genome>.txt: this is a parsed Phaster file, in csv format, containing columns contig, MGE and gene_nr.
"""

#import statements
import pandas as pd
from pathlib import Path
import csv
import gffutils
import sys
import os
import json

# script parameters
project_path=Path().resolve().parent.resolve().parent.resolve().parent
input_path=project_path / "results" / "intermediate" / "benchmarking" / "phaster_raw"
output_path=project_path / "results" / "intermediate" / "benchmarking" / "phaster_parsed"
gff_path=project_path / "data" / "legen_v4_dereplicated_gffs"

genome_file=sys.argv[1]
genome=str(genome_file.split('.')[0]+'.'+genome_file.split('.')[1])
print(genome)

def read_phaster(file, multi_contig=True):
    """Reads a phaster output file and generates a dataframe with this information."""
    df = pd.read_table(
    file,
    engine = 'c',
    skiprows =1,
    names = ['raw'],
    quoting=csv.QUOTE_NONE
    )
    
    # Select relevant rows only
    start=0
    for index,row in df.iterrows():
        if ("REGION" in str(row.item())):
            start=index 
    table_values=df.iloc[start:len(df)-1,:].reset_index(drop=True)
    
    # Separate into different columns (indicated by whitespaces)
    separated_table_values = pd.DataFrame.from_records(table_values.raw.apply(lambda s: s.split()))
    
    # First row contains the header
    separated_table_values.columns = separated_table_values.iloc[0,:]
    separated_table_values=separated_table_values.iloc[1:len(separated_table_values),:].reset_index(drop=True)
    
    if(not separated_table_values.empty):
        if(multi_contig):
            location=pd.DataFrame.from_records(separated_table_values.REGION_POSITION.apply(lambda s: s.split(sep=',')))
            separated_table_values['CONTIG']=location.iloc[:,0]
            separated_table_values['LOCATION_SEQUENCE']=" "
            for index,row in location.iterrows():
                for i in range(0,len(location.columns)):
                    if (not isinstance(location.iloc[index,i], type(None))):
                        if (":" in location.iloc[index,i]):
                            separated_table_values.LOCATION_SEQUENCE[index]=location.iloc[index,i]
            #separated_table_values['LOCATION_SEQUENCE']=str(location.iloc[:,-1])
            location2=pd.DataFrame.from_records(separated_table_values['LOCATION_SEQUENCE'].apply(lambda s: s.split(sep=':')))
            separated_table_values['REGION_POSITION']=location2.iloc[:,1]
        else:
            separated_table_values['CONTIG']= [1] * len(separated_table_values)


        # Add columns 'START' and 'END' (information extracted from 'REGION_POSITION')
        separated_table_values[['START', 'END']]=pd.DataFrame.from_records(separated_table_values.REGION_POSITION.apply(lambda s: s.split(sep='-')))

        return separated_table_values.loc[:, ['CONTIG','REGION', 'START', 'END']]
    else:
        return pd.DataFrame(columns=['CONTIG','REGION', 'START', 'END'])
    
def get_genes_multi_contig(df, start, stop, contig):
    """Gets the gene numbers corresponding to locations specified in basepairs, for a contig."""
    return df[(df.start>=start) & (df.end<=stop) & (df.contig==contig)]

def add_start_end_gene(genome, phaster, multi_contig=True):
    """Adds the columns 'START_GENE' and 'END_GENE' to the phaster file, indicating the location of prophages based on gene_nrs instead of basepairs."""
    if (not phaster.empty):
        gff= str(str(gff_path) + "/" + genome + ".gff")
        # generate a database
        gffutils.create_db(gff, str(genome + "_db"))
        db = gffutils.FeatureDB(dbfn=str(genome + "_db"))
        # generate a dataframe
        df=pd.DataFrame(columns = ['contig', 'ID', 'start', 'end', 'strand'])
        # fill up the dataframe, using the database
        query = db.execute("select seqid,start,end,strand,attributes from features where featuretype = 'CDS'")
        result = query.fetchall()

        # generate a dataframe
        gene_list=[]
        for each in result:
            gene_list.append({'ID':json.loads(each['attributes'])['ID'][0], 'contig':each['seqid'], 'end':each['end'], 'start':each['start'], 'strand':each['strand']})
        df=pd.DataFrame.from_records(gene_list)

        # add start genes
        start_genes=[]
        for index, row in phaster.iterrows():
            if (multi_contig):
                start_genes.append(get_genes_multi_contig(df, int(row.START), int(row.END), row.CONTIG).index[0])
            else:
                start_genes.append(get_genes(df, int(row.START), int(row.END)).index[0])
        phaster['START_GENE']=start_genes

        # add end genes
        end_genes=[]
        for index, row in phaster.iterrows():
            if (multi_contig):
                end_genes.append(get_genes_multi_contig(df, int(row.START), int(row.END), row.CONTIG).index[-1])
            else:
                end_genes.append(get_genes(df, int(row.START), int(row.END)).index[-1])
        phaster['END_GENE']=end_genes
        db=str(genome + "_db")
        os.remove(db)
        return phaster.loc[:,['CONTIG', 'REGION', 'START_GENE', 'END_GENE']]
    else:
        return pd.DataFrame(columns=['CONTIG', 'REGION', 'START_GENE', 'END_GENE'])

def transform_frame(phaster): 
    """Transforms the phaster dataframe, specifying the start and end genes of each prophage, into a dataframe where all genes belonging to a prophage are listed."""
    if (not phaster.empty):
        MGE_list=[]
        for index,row in phaster.iterrows():
            for gene in range(row.START_GENE, row.END_GENE+1):
                MGE_list.append({'contig':row.CONTIG, 'MGE':row.REGION, 'gene_nr':gene})
        MGE_frame=pd.DataFrame.from_records(MGE_list)
        return MGE_frame
    else:
        return pd.DataFrame(columns=['contig', 'MGE', 'gene_nr'])

def write_output(MGE_frame):
    """Writed the output into a csv file in the correct directory."""
    if (not MGE_frame.empty):
        MGE_frame.to_csv(output_path / genome_file)

phaster=read_phaster(input_path / genome_file)
phaster=add_start_end_gene(genome, phaster)
MGE_frame=transform_frame(phaster)
write_output(MGE_frame)
