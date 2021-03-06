# Thesis: Detection of novel mobile genetic elements in bacterial genomes

This code is meant for the detection of mobile genetic elements (MGEs) in bacterial genomes via two complementary comparative genomics strategies.

## Input data
The input data concerning a set of genomes of choice, consist of the following files:
- genomes_metadata.csv: this is a csv file containing the following columns: genome, gtdb_species, gtdb_genus, gtdb_family, quality, checkm_completeness, checkm_contamination, ncbi_isolation_source, ncbi_strain_identifiers, gtdb_representative,species. Only the first two columns are required.
- pangenome.tsv: this is a tsv file containing the following columns: gene_id, genome, orthogroup. Gene_id is a string that is built up as follows: <contig>.<rank>, where contig represents a unique identifier for a contig, and rank represents the order of the genes in this contig. This file can be generated by the SCARAP toolkit.
- \<genome\>.ffn: this is a nucleotide FASTA file for all genes in a genome. This data is required for all genomes.
- \<genome\>.gff: this is a general feature format for a genome. This data is only required for genomes for which the MGEs need to be detected.
    
## Code structure
The code is divided into the following directories:
- main code: most important code of the project, containing all steps to calculate phyletic distribution pattern based and phylogeny based mobility scores of orthorgroups and visualize them over genomes.
- training: code to train the different MGE detection strategies (fine-tuning of their parameters).
- benchmarking: code to measure performance of the different MGE detection strategies and code to parse the output of tools used for benchmarking (at present: CONJscan and Phaster).
- main notebooks: notebooks to generate plots and perform statistical analyses.
- extra code: temporary scripts and old scripts that have been replaced.
- exploratory notebooks: notebooks to try out code and for debugging.

