### Script to create barplots of taxonomy prevalence at various levels ###

# load pandas and matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# numpy is needed to count non zero, matplot.pyplot is needed for plot labels

# load kraken results table with all samples and all taxonomy levels
kraken_merged = pd.read_csv(snakemake.input[0], sep='\t')

kraken_kingdom = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                     .endswith('kraken')]].set_index(kraken_merged['kingdom']).drop('unidentified')
kraken_phylum = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                    .endswith('kraken')]].set_index(kraken_merged['phylum']).drop('unidentified')
kraken_class = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                   .endswith('kraken')]].set_index(kraken_merged['class']).drop('unidentified')
kraken_order = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                   .endswith('kraken')]].set_index(kraken_merged['order']).drop('unidentified')
kraken_family = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                    .endswith('kraken')]].set_index(kraken_merged['family']).drop('unidentified')
kraken_species = kraken_merged[kraken_merged.columns[pd.Series(kraken_merged.columns).str
                                                     .endswith('kraken')]].set_index(kraken_merged['species']).drop('unidentified')
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# barplot of top 25 most prevalent species from kraken results

kraken_species['prevalence'] = np.count_nonzero(kraken_species, axis=1)

# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_species_top25 = kraken_species.nlargest(25,'prevalence')
kraken_species_barplot = kraken_species_top25['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top 25 species prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Species")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[0])

# barplot of top 25 most prevalent families from kraken results

kraken_family_nonrepeat = kraken_family.reset_index().drop_duplicates(subset=['family']).set_index('family')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_family_nonrepeat['prevalence'] = np.count_nonzero(kraken_family_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_family_top25 = kraken_family_nonrepeat.nlargest(25,'prevalence')
kraken_family_barplot = kraken_family_top25['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top 25 family prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Families")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[1])

# barplot of top 25 most prevalent families from kraken results

kraken_class_nonrepeat = kraken_class.reset_index().drop_duplicates(subset=['class']).set_index('class')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_class_nonrepeat['prevalence'] = np.count_nonzero(kraken_class_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_class_top25 = kraken_class_nonrepeat.nlargest(25,'prevalence')
kraken_class_barplot = kraken_class_top25['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top 25 class prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Classes")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[2])

# barplot of top 25 most prevalent kingdoms from kraken results

kraken_kingdom_nonrepeat = kraken_kingdom.reset_index().drop_duplicates(subset=['kingdom']).set_index('kingdom')
# kraken results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

kraken_kingdom_nonrepeat['prevalence'] = np.count_nonzero(kraken_kingdom_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

kraken_kingdom_barplot = kraken_kingdom_nonrepeat['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top 25 kingdom prevalence from Kraken2 results")
plt.xlabel("Prevalence")
plt.ylabel("Kingdoms")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[3])
