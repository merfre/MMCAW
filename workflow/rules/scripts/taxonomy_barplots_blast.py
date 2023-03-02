### Script to create barplots of taxonomy prevalence at various levels ###

# load pandas and matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# numpy is needed to count non zero, matplot.pyplot is needed for plot labels

prevno = snakemake.config["prevalence"]

# load blast results table with all samples and all taxonomy levels
blast_merged = pd.read_csv(snakemake.input[0], sep='\t')

blast_kingdom = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                     .endswith('blast')]].set_index(blast_merged['kingdom']).drop('u')
blast_phylum = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                    .endswith('blast')]].set_index(blast_merged['phylum']).drop('unidentified')
blast_class = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                   .endswith('blast')]].set_index(blast_merged['class']).drop('unidentified')
blast_order = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                   .endswith('blast')]].set_index(blast_merged['order']).drop('unidentified')
blast_family = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                    .endswith('blast')]].set_index(blast_merged['family']).drop('unidentified')
blast_species = blast_merged[blast_merged.columns[pd.Series(blast_merged.columns).str
                                                     .endswith('blast')]].set_index(blast_merged['species']).drop('unidentified')
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# barplot of top most prevalent species from blast results

blast_species['prevalence'] = np.count_nonzero(blast_species, axis=1)

# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

blast_species_top = blast_species.nlargest(prevno,'prevalence')
blast_species_barplot = blast_species_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title(print("Barplot of top " + str(prevno) + " species prevalence from BLAST results"))
plt.xlabel("Prevalence")
plt.ylabel("Species")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[0])

# barplot of top most prevalent families from blast results

blast_family_nonrepeat = blast_family.reset_index().drop_duplicates(subset=['family']).set_index('family')
# blast results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

blast_family_nonrepeat['prevalence'] = np.count_nonzero(blast_family_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

blast_family_top = blast_family_nonrepeat.nlargest(prevno,'prevalence')
blast_family_barplot = blast_family_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title(print("Barplot of top " + str(prevno) + " family prevalence from BLAST results"))
plt.xlabel("Prevalence")
plt.ylabel("Families")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[1])

# barplot of top most prevalent families from blast results

blast_class_nonrepeat = blast_class.reset_index().drop_duplicates(subset=['class']).set_index('class')
# blast results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

blast_class_nonrepeat['prevalence'] = np.count_nonzero(blast_class_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

blast_class_top = blast_class_nonrepeat.nlargest(prevno,'prevalence')
blast_class_barplot = blast_class_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title(print("Barplot of top " + str(prevno) + " class prevalence from BLAST results"))
plt.xlabel("Prevalence")
plt.ylabel("Classes")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[2])

# barplot of top most prevalent kingdoms from blast results

blast_kingdom_nonrepeat = blast_kingdom.reset_index().drop_duplicates(subset=['kingdom']).set_index('kingdom')
# blast results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

blast_kingdom_nonrepeat['prevalence'] = np.count_nonzero(blast_kingdom_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

blast_kingdom_barplot = blast_kingdom_nonrepeat['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title(print("Barplot of top " + str(prevno) + " kingdom prevalence from BLAST results"))
plt.xlabel("Prevalence")
plt.ylabel("Kingdoms")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[3])
