### Script to create barplots of taxonomy prevalence at various levels ###

# load pandas and matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# numpy is needed to count non zero, matplot.pyplot is needed for plot labels

prevno = snakemake.config["prevalence"]

# load cat results table with all samples and all taxonomy levels
cat_merged = pd.read_csv(snakemake.input[0], sep='\t')

cat_kingdom = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                     .endswith('cat')]].set_index(cat_merged['kingdom']).drop('unidentified')
cat_phylum = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                    .endswith('cat')]].set_index(cat_merged['phylum']).drop('unidentified')
cat_class = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                   .endswith('cat')]].set_index(cat_merged['class']).drop('unidentified')
cat_order = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                   .endswith('cat')]].set_index(cat_merged['order']).drop('unidentified')
cat_family = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                    .endswith('cat')]].set_index(cat_merged['family']).drop('unidentified')
cat_species = cat_merged[cat_merged.columns[pd.Series(cat_merged.columns).str
                                                     .endswith('cat')]].set_index(cat_merged['species']).drop('unidentified')
# to remove unclassified reads rows with u, unknown, or unidentified unidentified assignment

# barplot of top most prevalent species from cat results

cat_species['prevalence'] = np.count_nonzero(cat_species, axis=1)

# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

cat_species_top = cat_species.nlargest(prevno,'prevalence')
cat_species_barplot = cat_species_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " species prevalence from CAT results")
plt.xlabel("Prevalence")
plt.ylabel("Species")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[0])

# barplot of top most prevalent families from cat results

cat_family_nonrepeat = cat_family.reset_index().drop_duplicates(subset=['family']).set_index('family')
# cat results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

cat_family_nonrepeat['prevalence'] = np.count_nonzero(cat_family_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

cat_family_top = cat_family_nonrepeat.nlargest(prevno,'prevalence')
cat_family_barplot = cat_family_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " family prevalence from CAT results")
plt.xlabel("Prevalence")
plt.ylabel("Families")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[1])

# barplot of top most prevalent families from cat results

cat_class_nonrepeat = cat_class.reset_index().drop_duplicates(subset=['class']).set_index('class')
# cat results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

cat_class_nonrepeat['prevalence'] = np.count_nonzero(cat_class_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

cat_class_top = cat_class_nonrepeat.nlargest(prevno,'prevalence')
cat_class_barplot = cat_class_top['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " class prevalence from CAT results")
plt.xlabel("Prevalence")
plt.ylabel("Classes")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[2])

# barplot of top most prevalent kingdoms from cat results

cat_kingdom_nonrepeat = cat_kingdom.reset_index().drop_duplicates(subset=['kingdom']).set_index('kingdom')
# cat results list the taxonomies multiple times and duplicates (in the index) must be removed to assess prevalence

cat_kingdom_nonrepeat['prevalence'] = np.count_nonzero(cat_kingdom_nonrepeat, axis=1)
# count nonzero of rows will determine the prevalence of the taxonomy and that's added as a column to the dataframe

cat_kingdom_barplot = cat_kingdom_nonrepeat['prevalence'].plot(kind="barh")
# change color with: color=[''], barh makes plot horizontal
plt.title("Barplot of top " + prevno + " kingdom prevalence from CAT results")
plt.xlabel("Prevalence")
plt.ylabel("Kingdoms")
plt.tight_layout()
# tight layout needed to not cut off species names
plt.savefig(snakemake.output[3])
