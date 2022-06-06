#!/usr/bin/env python
# coding: utf-8

"""Visualizes and records mobility scores of orthogroups over a genome.

Visualizes and records the mobility score of orthogroups over a genome of choice to detect mobile genetic elements (MGEs).

Input files: 
genomes_metadata.csv: this is a csv file containing the following columns: genome, gtdb_species, gtdb_genus, gtdb_family, quality, checkm_completeness, checkm_contamination, ncbi_isolation_source, ncbi_strain_identifiers, gtdb_representative,species. Only the first two columns are required.
mobility_frame.csv: this is a csv file containing the following columns: orthogroup, species_count, genomes_count, accessory_fraction.
pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>_<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig.

Output files: 
<genome>-<contig>.png: mobility plot showing the mobility scores for cosecutive genes of a genome of choice.
<genome>-<contig>.csv: mobility frame for consecutive genes of a genome of choice.
"""

# import statements
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import sys
from matplotlib.lines import Line2D
import os

output_dir=sys.argv[2]

# script parameters
genomes_per_species = 10 # minimal number of genomes per species
threshold_accessory = 0.9 # fraction of genomes per species an orthogroup can be present in at max. to be considered accessory
min_contigsize = 20 # minimal contig size to be able to detect MGEs
small_contig_threshold = 100 # contigs containing a number of genes >= 100 are considered small
running_average_default = 20 # default number of consecutive genes to take into account to calculate running averages
running_average_default_small_contigs = 5 # default number of consecutive genes to take into account to calculate running averages in small contigs
project_path = Path().resolve().parent.parent
path_genomes = project_path / "data" / "genomes_metadata.csv"
path_mobility = project_path / "results" / "intermediate" / "mobility_frame.csv"
path_pangenome = project_path / "data" / "pangenome.tsv"
output_path_plots = project_path / "results" / "mobility_plots" / output_dir
output_path_plots_2 = project_path / "results" / "tree_mobility_plots"
output_path_files = project_path / "results" / "mobility_files" / output_dir

if (os.path.exists(output_path_plots) == False):
    os.mkdir(output_path_plots)
if (os.path.exists(output_path_files) == False):
    os.mkdir(output_path_files)

# import data
genomes = pd.read_csv(path_genomes).loc[:,'genome':'gtdb_species']
mobility_frame = pd.read_csv(path_mobility)
pangenome = pd.read_csv(path_pangenome, delimiter="\t", header=None)
pangenome.columns = ["gene", "genome", "orthogroup"]

def get_species(genomes, genome):
    """Returns the bacterial species to which a genome belongs, based on a genomes metadata file."""
    return genomes[genomes.genome==genome].gtdb_species.item()

def determine_order(pangenome):
    """Extracts contig and order information from the pangenome file."""
    order = pangenome["gene"].str.split("_", n=1, expand=True)
    order[1] = pd.to_numeric(order[1])
    order.columns = ["contig", "order"]
    full_order = order.merge(pangenome, left_index=True, right_index=True)
    return full_order

def get_nr_genes(full_order):
    """Counts the number of genes present per contig in a dataframe with columns 'contig' and 'gene'."""
    nr_genes = full_order.reindex(columns=["contig", "gene"]).drop_duplicates().groupby(by = "contig", as_index = False).count()
    nr_genes.columns = ["contig", "nr_genes"]
    return nr_genes

def get_size(full_order, contig):
    """Calculates the size of a contig of choice based on a dataframe with columns 'contig' and 'gene'."""
    nr_genes = get_nr_genes(full_order)
    return nr_genes.nr_genes[nr_genes.contig==contig].item()

def get_contigs(full_order, genome, filter_size=min_contigsize):
    """Returns the contigs belonging to a genome of choice, in a dataframe with columns 'contig' and 'gene', filtered to contain a minimal number of genes per contig."""
    contigs = pd.DataFrame()
    contigs['contig']=full_order.contig[full_order.genome==genome].drop_duplicates()
    contigs_size=contigs.merge(get_nr_genes(full_order), on='contig')
    contigs_final=contigs_size[(contigs_size.nr_genes)>=filter_size] # filter contigs: only contigs exceeding a size limit are retained
    return contigs_final.contig

def count_genomes_per_species(genomes):
    """Counts the number of genomes present per species present in a genomes metadata file."""
    counts = genomes.groupby(by="gtdb_species", as_index=False).count()
    counts.columns = ['gtdb_species', 'counts']
    return counts

def filter_species(genomes, filter_threshold=genomes_per_species):
    """Filters the species to contain a minimal number of genomes per species."""
    counts = count_genomes_per_species(genomes)
    genomes_filtered = pd.merge(genomes, counts[counts.counts >= filter_threshold], on='gtdb_species', how='right').sort_values('gtdb_species')
    return genomes_filtered


def get_scc(genomes, genome, full):
    """Calculates flags for orthogroups present in a genome of choice:
    whether or not this orthogroup is present as an accessory gene (represented by accessory==True), and whether or not this orthogroup is present as a single-copy gene (reprsented by 'count'==1). Returns a dataframe with columns 'orthogroup', 'count', 'accessory' and 'scc'. 'scc' is True when the orthogroup is accessory in this genome or when it is present in multiple copies.
    
    Args:
        genomes: Genomes metadata file
        genome: String representing the ID of the genome for which the mobility score needs to be plotted.
        full: Datafame containing pangenome information, filtered for relevant species.
    """
    
    species = get_species(genomes, genome)
    counts = count_genomes_per_species(genomes)
    full_species = full[full.gtdb_species==species]
    orthocounts_genomes_per_species = full_species.loc[:,["orthogroup","gtdb_species","genome"]].drop_duplicates().groupby(by = ["gtdb_species", "orthogroup"], as_index = False).count()
    orthocounts_genomes_per_species.columns = ["gtdb_species", "orthogroup", "orthocounts"]
    res = pd.merge(orthocounts_genomes_per_species, counts, on = 'gtdb_species', how = 'left')
    res.columns = ["gtdb_species", "orthogroup", "orthocounts", "genome_counts"]
    res['accessory'] = (res.orthocounts <= threshold_accessory * res.genome_counts)
    
    full_species = full[full.gtdb_species==species]
    orthocounts_genomes_per_species = full_species.loc[:,["orthogroup","gtdb_species","genome"]].drop_duplicates().groupby(by = ["gtdb_species", "orthogroup"], as_index = False).count()
    orthocounts_genomes_per_species.columns = ["gtdb_species", "orthogroup", "orthocounts"]
    res = pd.merge(orthocounts_genomes_per_species, counts, on = 'gtdb_species', how = 'left')
    res.columns = ["gtdb_species", "orthogroup", "orthocounts", "genome_counts"]
    res['accessory'] = (res.orthocounts <= threshold_accessory * res.genome_counts)
    
    full_genome = full.loc[full.genome==genome]
    genome_orthocounts = full_genome.loc[:,["orthogroup","gtdb_species"]].groupby(by = ["orthogroup"], as_index = False).count()
    genome_orthocounts.columns = ["orthogroup", "count"]
    
    core_copy = pd.merge(genome_orthocounts.loc[:,["orthogroup", "count"]], res.loc[:,["orthogroup", "accessory"]], on="orthogroup", how="left")
    d = {True: 'red', False: 'blue'}
    core_copy['scc']=(core_copy.accessory)|(core_copy['count']>1)
    core_copy.scc = core_copy.scc.map(d)
    
    return core_copy

def mobility_plots(genomes, genome, full, order, method="mobility_score", n=running_average_default, indicate_accessory=True):
    """Plots the mobility score across a genome of choice, according to a method of choice.
    
    Args:
        genomes: Genomes metadata file
        genome: String representing the ID of the genome for which the mobility score needs to be plotted.
        full: Datafame containing pangenome information, filtered for relevant species.
        order: Datafame containing information on the contig, and order of genes, and pangenome information.
        method: String representing the method used for calculating the mobility score. Options: "species_count", "genomes_count", "accessory_fraction".
        n: integer indicating the number of consecutive genes to take into account for a moving average calulation.
        indicate_accessory: Boolean indicating whether or not the plots need to be coloured based on accessory and copy number information of orthogroups.
    """
    for contig in get_contigs(order, genome):
        # a different number of consecutive genes is taken into account for small contigs
        if get_size(full_order, contig) <= small_contig_threshold:
            n=running_average_default_small_contigs
        scc=get_scc(genomes, genome, full)
        genome_order = order[(order.genome == genome)&(order.contig ==contig)].loc[:,["order", "orthogroup"]]
        mob_order = pd.merge(genome_order, mobility_frame, on="orthogroup", how="left").sort_values(by="order")
        mob_order_scc = pd.merge(mob_order, scc, on="orthogroup", how = "left")

        # plot 1
        fig, ax = plt.subplots(figsize=(30, 10))
        y=mob_order_scc[method].rolling(window=n, center=True).mean()
        ax.plot(mob_order_scc.order, y)
        # color points based on acesory and copy number information of orthogroups
        #if indicate_accessory:
        #    ax.scatter(mob_order_scc.order, mob_order_scc[method].rolling(window=n, center=True).mean(), c=mob_order_scc.scc, marker='o')
        if (contig=="MLNO01000060.1"):
            x=list(range(71,84))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MLNO01000026.1"):
            x=list(range(0,23))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MLNO01000008.1"):
            x=list(range(0,27))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MLNO01000011.1"):
            x=list(range(10,20))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MLNO01000070.1"):
            x=list(range(2, 12))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MULL01000001.1"):
            #x=list(range(707, 723))
            #x.extend(range(731,751))
            #x.extend(range(1149,1160))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster predictions")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster predictions')]
            #ax.legend(handles=legend_elements, fontsize=20)
            
            x=list(range(1,36))
            x.extend(range(702, 725))
            x.extend(range(732, 753))
            x.extend(range(2255, 2300))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE predictions by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE predictions by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

        elif (contig=="MULL01000002.1"):
            #x=list(range(11, 34))
            #x.extend(range(27,49))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)
            
            x=list(range(0,56))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(36, 57))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the detailed approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the detailed approach')]
            #ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000003.1"):
            #x=list(range(20, 30))
            #x.extend(range(46,65))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster predictions")

            #x2=[1, 2, 3, 7, 38]
            #y2=y.iloc[x2]
            #ax.scatter([i+1 for i in x2], y2, c='blue', marker='o', label = "CONJscan prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster predictions'), Line2D([0], [0], marker='o', color='blue', label='CONJscan prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)


            x=list(range(0,66))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(1, 26))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the simple approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the simple approach')]
            #ax.legend(handles=legend_elements, fontsize=20)


            #x=list(range(1, 26))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the detailed approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the detailed approach')]
            #ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000004.1"):
            #x=list(range(4, 16))
            #x.extend(range(34, 47))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster predictions")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster predictions')]
            #ax.legend(handles=legend_elements, fontsize=20)

            x=list(range(0,47))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000005.1"):
            #x=list(range(20, 29))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)

            x=list(range(0,48))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000006.1"):
            #x=list(range(0, 17))
            #x.extend(range(35, 51))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster predictions")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster predictions')]
            #ax.legend(handles=legend_elements, fontsize=20)

            x=list(range(0,53))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(30,51))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the simple approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the simple approach')]
            #ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(30, 51))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the detailed approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the detailed approach')]
            #ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000007.1"):
            #x=list(range(0, 23))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)

            x=list(range(0,38))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000009.1"):
            #x=[17, 21, 22, 24, 28, 29, 31, 5, 8, 9, 10, 11, 15, 16, 9]
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "CONJscan prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='CONJscan prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)

            x=list(range(0,87))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(2,20))
            #x.extend(range(35, 64))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the simple approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the simple approach')]
            #ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(1, 30))
            #x.extend(range(36, 79))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the detailed approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the detailed approach')]
            #ax.legend(handles=legend_elements, fontsize=20)


        elif (contig=="MULL01000010.1"):
            x=list(range(0,40))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

            #x=list(range(6, 42))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "MGE prediction by the detailed approach")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='MGE prediction by the detailed approach')]
            #ax.legend(handles=legend_elements, fontsize=20)



        #elif (contig=="MPAQ01000021.1"):
        #    x=list(range(2, 18))
        #    y1=y.iloc[x]
        #    ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
        #    legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
        #    ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MPAQ01000049.1"):
            x=list(range(4, 13))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)
        elif (contig=="MPAQ01000005.1"):
            #x=list(range(6, 69))
            #y1=y.iloc[x]
            #ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            #legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            #ax.legend(handles=legend_elements, fontsize=20)
            
            x=list(range(22, 72))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='green', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='green', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)

        elif (contig=="MPAQ01000079.1"):
            x=list(range(0, 7))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='red', marker='o', label = "Phaster prediction")
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)

        elif (contig=="MPAQ01000021.1"):
            x=list(range(0, 24))
            y1=y.iloc[x]
            ax.scatter([i+1 for i in x], y1, c='green', marker='o', label = "MGE prediction by the baseline approach")
            legend_elements = [Line2D([0], [0], marker='o', color='green', label='MGE prediction by the baseline approach')]
            ax.legend(handles=legend_elements, fontsize=20)









        ax.set_xlabel("Gene position", fontsize=20)
        ax.set_ylabel("Mobility score", fontsize=20)
        ax.set_title(genome, fontsize=30)
        ax.tick_params(labelsize=18)

        if ((contig=="MLNO01000060.1")|(contig=="MLNO01000011.1")):
            legend_elements = [Line2D([0], [0], marker='o', color='red', label='Phaster prediction')]
            ax.legend(handles=legend_elements, fontsize=20)

        #legend_elements = [Line2D([0], [0], marker='o', color='blue', label='Single copy core'),
        #           Line2D([0], [0], marker='o', color='red', label='Accessory and/or multi copy')]
                #plt.savefig(output_path_plots / str(genome+'-'+contig+'.png'))
        plt.savefig(output_path_plots / str("baseline_"+genome+'-'+contig+'.png'))


        # plot 2
        fig, ax = plt.subplots(figsize=(30, 10))
        ax.plot(mob_order_scc.order, mob_order_scc['tree_score'].rolling(window=n, center=True).mean())
        # color points based on acesory and copy number information of orthogroups
        if indicate_accessory:
            ax.scatter(mob_order_scc.order, mob_order_scc['tree_score'].rolling(window=n, center=True).mean(), c=mob_order_scc.scc, marker='o')
        ax.set_xlabel("Gene position", fontsize=20)
        ax.set_ylabel("Mobility score", fontsize=20)
        ax.set_title(genome, fontsize=30)
        ax.tick_params(labelsize=18)
        legend_elements = [Line2D([0], [0], marker='o', color='blue', label='Single copy core'),
                   Line2D([0], [0], marker='o', color='red', label='Accessory and/or multi copy')]
        ax.legend(handles=legend_elements, fontsize=20)
        plt.savefig(output_path_plots_2 / str(genome+'-'+contig+'.png'))
        
def mobility_file(genomes, genome, full, order):
    """Records the mobility score across a genome of choice in a csv file.
    
    Args:
        genomes: Genomes metadata file
        genome: String representing the ID of the genome for which the mobility score needs to be plotted.
        full: Datafame containing pangenome information, filtered for relevant species.
        order: Datafame containing information on the contig, and order of genes, and pangenome information.
    """
    for contig in get_contigs(order, genome):
        scc=get_scc(genomes, genome, full)
        genome_order = order[(order.genome == genome)&(order.contig ==contig)].loc[:,["order", "orthogroup", "gene"]]
        mob_order = pd.merge(genome_order, mobility_frame.loc[:,["orthogroup", "accessory_fraction", "tree_score"]], on="orthogroup", how="left").sort_values(by="order")
        mob_order_scc = pd.merge(mob_order, scc.loc[:, ["orthogroup", "count", "accessory"]], on="orthogroup", how = "left")
        mob_order_scc.to_csv(output_path_files / str(genome+'-'+contig+'.csv'))
    
genomes_filtered = filter_species(genomes=genomes)
full = pd.merge(genomes_filtered, pangenome, on='genome', how='left')
full_order = determine_order(pangenome)
mobility_plots(genomes=genomes, genome=sys.argv[1], full=full, order=full_order, method="accessory_fraction")
mobility_file(genomes=genomes, genome=sys.argv[1], full=full, order=full_order)
