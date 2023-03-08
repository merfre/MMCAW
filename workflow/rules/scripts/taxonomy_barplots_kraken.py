### Script to create barplots of taxonomy prevalence at various levels ###

# load pandas and matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# numpy is needed to count non zero, matplot.pyplot is needed for plot labels

prevno = snakemake.config["prevalence"]

# load kraken results table with all samples and all taxonomy levels
kraken_merged = pd.read_csv(snakemake.input[0], sep='\t')

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

# barplot of top most prevalent species from kraken results

kraken_species['prevalence'] = np.count_nonzero(kraken_species, axis=1)

# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_species_top = kraken_species.nlargest(int(prevno),'prevalence')
kraken_species_barplot = kraken_species_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " species prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Species")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[0])

# barplot of top most prevalent families from kraken results

kraken_family_nonrepeat = kraken_family.reset_index().drop_duplicates(subset=['family']).set_index('family')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_family_nonrepeat['prevalence'] = np.count_nonzero(kraken_family_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_family_top = kraken_family_nonrepeat.nlargest(int(prevno),'prevalence')
kraken_family_barplot = kraken_family_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " family prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Families")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[1])

# barplot of top most prevalent families from kraken results

kraken_class_nonrepeat = kraken_class.reset_index().drop_duplicates(subset=['class']).set_index('class')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_class_nonrepeat['prevalence'] = np.count_nonzero(kraken_class_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_class_top = kraken_class_nonrepeat.nlargest(int(prevno),'prevalence')
kraken_class_barplot = kraken_class_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " class prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Classes")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[2])

# barplot of top most prevalent kingdoms from kraken results

kraken_kingdom_nonrepeat = kraken_kingdom.reset_index().drop_duplicates(subset=['kingdom']).set_index('kingdom')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_kingdom_nonrepeat['prevalence'] = np.count_nonzero(kraken_kingdom_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_kingdom_barplot = kraken_kingdom_nonrepeat['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " kingdom prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Kingdoms")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[3])
