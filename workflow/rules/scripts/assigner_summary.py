### Script for creating a summary table of taxonomy assigner results and performance ###

import pandas as pd
import numpy as np
# numpy is needed to count non zero

# load cat results table with all samples and all taxonomy levels
cat_merged = pd.read_csv(snakemake.input[0], sep='\t')

# separate results at every taxonomy level with taxon as rows and samples as columns
cat_kingdom = cat_merged[~cat_merged.kingdom.str.fullmatch("unidentified")]
cat_kingdom= cat_kingdom[cat_kingdom.columns[pd.Series(cat_kingdom.columns).str
                                                     .endswith('cat')]].set_index(cat_kingdom['kingdom'])

cat_phylum = cat_merged[~cat_merged.phylum.str.fullmatch("unidentified")]
cat_phylum = cat_phylum[cat_phylum.columns[pd.Series(cat_phylum.columns).str
                                                    .endswith('cat')]].set_index(cat_phylum['phylum'])

cat_class = cat_merged[~cat_merged['class'].str.fullmatch("unidentified")]
cat_class = cat_class[cat_class.columns[pd.Series(cat_class.columns).str
                                                   .endswith('cat')]].set_index(cat_class['class'])

cat_order = cat_merged[~cat_merged.order.str.fullmatch("unidentified")]
cat_order = cat_order[cat_order.columns[pd.Series(cat_order.columns).str
                                                   .endswith('cat')]].set_index(cat_order['order'])

cat_family = cat_merged[~cat_merged.family.str.fullmatch("unidentified")]
cat_family = cat_family[cat_family.columns[pd.Series(cat_family.columns).str
                                                    .endswith('cat')]].set_index(cat_family['family'])

cat_species = cat_merged[~cat_merged.species.str.fullmatch("unidentified")]
cat_species= cat_species[cat_species.columns[pd.Series(cat_species.columns).str
                                                     .endswith('cat')]].set_index(cat_species['species'])
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# load kraken2 results table with all samples and all taxonomy levels
kraken_merged = pd.read_csv(snakemake.input[1], sep='\t')

kraken_kingdom = kraken_merged[~kraken_merged.kingdom.str.fullmatch("unidentified")]
kraken_kingdom= kraken_kingdom[kraken_kingdom.columns[pd.Series(kraken_kingdom.columns).str
                                                     .endswith('kraken')]].set_index(kraken_kingdom['kingdom'])

kraken_phylum = kraken_merged[~kraken_merged.phylum.str.fullmatch("unidentified")]
kraken_phylum = kraken_phylum[kraken_phylum.columns[pd.Series(kraken_phylum.columns).str
                                                    .endswith('kraken')]].set_index(kraken_phylum['phylum'])

kraken_class = kraken_merged[~kraken_merged['class'].str.fullmatch("unidentified")]
kraken_class = kraken_class[kraken_class.columns[pd.Series(kraken_class.columns).str
                                                   .endswith('kraken')]].set_index(kraken_class['class'])

kraken_order = kraken_merged[~kraken_merged.order.str.fullmatch("unidentified")]
kraken_order = kraken_order[kraken_order.columns[pd.Series(kraken_order.columns).str
                                                   .endswith('kraken')]].set_index(kraken_order['order'])

kraken_family = kraken_merged[~kraken_merged.family.str.fullmatch("unidentified")]
kraken_family = kraken_family[kraken_family.columns[pd.Series(kraken_family.columns).str
                                                    .endswith('kraken')]].set_index(kraken_family['family'])

kraken_species = kraken_merged[~kraken_merged.species.str.fullmatch("unidentified")]
kraken_species= kraken_species[kraken_species.columns[pd.Series(kraken_species.columns).str
                                                     .endswith('kraken')]].set_index(kraken_species['species'])
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# load blast results table with all samples and all taxonomy levels
blast_merged = pd.read_csv(snakemake.input[2], sep='\t')

blast_kingdom = blast_merged[~blast_merged.kingdom.str.fullmatch("u")]
blast_kingdom= blast_kingdom[blast_kingdom.columns[pd.Series(blast_kingdom.columns).str
                                                     .endswith('blast')]].set_index(blast_kingdom['kingdom'])

blast_phylum = blast_merged[~blast_merged.phylum.str.fullmatch("unidentified")]
blast_phylum = blast_phylum[blast_phylum.columns[pd.Series(blast_phylum.columns).str
                                                    .endswith('blast')]].set_index(blast_phylum['phylum'])

blast_class = blast_merged[~blast_merged['class'].str.fullmatch("unidentified")]
blast_class = blast_class[blast_class.columns[pd.Series(blast_class.columns).str
                                                   .endswith('blast')]].set_index(blast_class['class'])

blast_order = blast_merged[~blast_merged.order.str.fullmatch("unidentified")]
blast_order = blast_order[blast_order.columns[pd.Series(blast_order.columns).str
                                                   .endswith('blast')]].set_index(blast_order['order'])

blast_family = blast_merged[~blast_merged.family.str.fullmatch("unidentified")]
blast_family = blast_family[blast_family.columns[pd.Series(blast_family.columns).str
                                                    .endswith('blast')]].set_index(blast_family['family'])

blast_species = blast_merged[~blast_merged.species.str.fullmatch("unidentified")]
blast_species= blast_species[blast_species.columns[pd.Series(blast_species.columns).str
                                                     .endswith('blast')]].set_index(blast_species['species'])
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# create a taxonomy assigner summary table

Taxonomy_assigner_summary = pd.DataFrame(index=["Average reads assigned to species","Average percent of reads assigned to species","Average species richness","Most common species identified",
                                               "Average reads assigned to family","Average family richness","Most common family identified",
                                               "Average reads assigned to class","Average class richness","Most common class identified",
                                               "Average reads assigned to kingdom","Average kingdom richness","Most common kingdom identified"],
                                         columns=["CAT","Kraken2","BLAST"])

# number of reads assigned to the species level, on average

Taxonomy_assigner_summary.at['Average reads assigned to species', 'CAT'] = cat_species.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to species', 'Kraken2'] = kraken_species.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to species', 'BLAST'] = blast_species.sum().mean()
# sum of reads in each column then mean of all columns

# percent of reads assigned to the species level (out of total number of reads put into assigner), on average

seq_stats = pd.read_csv(snakemake.params[0], sep='\t')
mean_sequences = seq_stats['num_seqs'].mean()
# load the statistics calculated by fastqc then get the average number of sequences in each sample that went into the assigners

Taxonomy_assigner_summary.at['Average percent of reads assigned to species', 'CAT'] = cat_species.sum().mean() / mean_sequences * 100
Taxonomy_assigner_summary.at['Average percent of reads assigned to species', 'Kraken2'] = kraken_species.sum().mean() / mean_sequences * 100
Taxonomy_assigner_summary.at['Average percent of reads assigned to species', 'BLAST'] = blast_species.sum().mean() / mean_sequences * 100
# divide the value calculated previously by the average number of sequences and multiple by 100 for percentage

# number of different species identified, on average - average species richness

Taxonomy_assigner_summary.at['Average species richness', 'CAT'] = np.count_nonzero(cat_species, axis=0).mean()
Taxonomy_assigner_summary.at['Average species richness', 'Kraken2'] = np.count_nonzero(kraken_species.reset_index().drop_duplicates(subset=['species']), axis=0).mean()
Taxonomy_assigner_summary.at['Average species richness', 'BLAST'] = np.count_nonzero(blast_species.reset_index().drop_duplicates(subset=['species']), axis=0).mean()
# count the number of species with at least one read for each column then get mean of all columns
# kraken2 results list the taxonomies multiple times and duplicates (in the index) must be removed to assess sample richness

# most common species reads are assigned to

Taxonomy_assigner_summary.at['Most common species identified', 'CAT'] = cat_species.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common species identified', 'Kraken2'] = kraken_species.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common species identified', 'BLAST'] = blast_species.sum(axis=1).idxmax()
# get the sum of each row then get the index of the rox with highest number of reads

# number of reads assigned to the family level, on average

Taxonomy_assigner_summary.at['Average reads assigned to family', 'CAT'] = cat_family.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to family', 'Kraken2'] = kraken_family.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to family', 'BLAST'] = blast_family.sum().mean()
# sum of reads in each column then mean of all columns

# number of different families identified, on average

Taxonomy_assigner_summary.at['Average family richness', 'CAT'] = np.count_nonzero(cat_family, axis=0).mean()
Taxonomy_assigner_summary.at['Average family richness', 'Kraken2'] = np.count_nonzero(kraken_family.reset_index().drop_duplicates(subset=['family']), axis=0).mean()
Taxonomy_assigner_summary.at['Average family richness', 'BLAST'] = np.count_nonzero(blast_family.reset_index().drop_duplicates(subset=['family']), axis=0).mean()
# count the number of species with at least one read for each column then get mean of all columns
# kraken2 results list the taxonomies multiple times and duplicates (in the index) must be removed to assess sample richness

# most common family reads are assigned to

Taxonomy_assigner_summary.at['Most common family identified', 'CAT'] = cat_family.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common family identified', 'Kraken2'] = kraken_family.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common family identified', 'BLAST'] = blast_family.sum(axis=1).idxmax()
# get the sum of each row then get the index of the rox with highest number of reads

# number of reads assigned to the class level, on average

Taxonomy_assigner_summary.at['Average reads assigned to class', 'CAT'] = cat_class.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to class', 'Kraken2'] = kraken_class.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to class', 'BLAST'] = blast_class.sum().mean()
# sum of reads in each column then mean of all columns

# number of different classes identified, on average

Taxonomy_assigner_summary.at['Average class richness', 'CAT'] = np.count_nonzero(cat_class, axis=0).mean()
Taxonomy_assigner_summary.at['Average class richness', 'Kraken2'] = np.count_nonzero(kraken_class.reset_index().drop_duplicates(subset=['class']), axis=0).mean()
Taxonomy_assigner_summary.at['Average class richness', 'BLAST'] = np.count_nonzero(blast_class.reset_index().drop_duplicates(subset=['class']), axis=0).mean()
# count the number of species with at least one read for each column then get mean of all columns
# kraken2 results list the taxonomies multiple times and duplicates (in the index) must be removed to assess sample richness

# most common class reads are assigned to

Taxonomy_assigner_summary.at['Most common class identified', 'CAT'] = cat_class.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common class identified', 'Kraken2'] = kraken_class.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common class identified', 'BLAST'] = blast_class.sum(axis=1).idxmax()
# get the sum of each row then get the index of the rox with highest number of reads

# number of reads assigned to the kingdom level, on average

Taxonomy_assigner_summary.at['Average reads assigned to kingdom', 'CAT'] = cat_kingdom.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to kingdom', 'Kraken2'] = kraken_kingdom.sum().mean()
Taxonomy_assigner_summary.at['Average reads assigned to kingdom', 'BLAST'] = blast_kingdom.sum().mean()
# sum of reads in each column then mean of all columns

# number of different kingdoms identified, on average

Taxonomy_assigner_summary.at['Average kingdom richness', 'CAT'] = np.count_nonzero(cat_kingdom, axis=0).mean()
Taxonomy_assigner_summary.at['Average kingdom richness', 'Kraken2'] = np.count_nonzero(kraken_kingdom.reset_index().drop_duplicates(subset=['kingdom']), axis=0).mean()
Taxonomy_assigner_summary.at['Average kingdom richness', 'BLAST'] = np.count_nonzero(blast_kingdom.reset_index().drop_duplicates(subset=['kingdom']), axis=0).mean()
# count the number of species with at least one read for each column then get mean of all columns
# kraken2 results list the taxonomies multiple times and duplicates (in the index) must be removed to assess sample richness

# most common kindom reads are assigned to

Taxonomy_assigner_summary.at['Most common kingdom identified', 'CAT'] = cat_kingdom.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common kingdom identified', 'Kraken2'] = kraken_kingdom.sum(axis=1).idxmax()
Taxonomy_assigner_summary.at['Most common kingdom identified', 'BLAST'] = blast_kingdom.sum(axis=1).idxmax()
# get the sum of each row then get the index of the rox with highest number of reads

Taxonomy_assigner_summary.to_csv(snakemake.output[0],sep="\t")
