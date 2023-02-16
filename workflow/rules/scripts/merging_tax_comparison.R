### Script to combine and merge tax comparison summary results into 1 table ###

library(dplyr)
# Package required to get distinct taxonomies and levels
library(tidyr)
# Package to drop NAs

file_paths <- t(read.csv(snakemake@input[[1]], header=FALSE))
# Load a list of file paths of the individual CAT results

merge_table <- data.frame(matrix(ncol = 11, nrow = 1))
# create an empty dataframe to merge with the sample summaries of taxonomy comparison

colnames(merge_table)<- c("sample_id", "kingdom_no_agreements",	"phylum_no_agreements",	"class_no_agreements",	"order_no_agreements",	"family_no_agreements",	"genus_no_agreements",	"species_no_agreements",	"blast_cat_agreement",	"kraken_blast_agreement",	"cat_kraken_agreement")
# name the columns the same as the column names of the summary tables from tax comparison

for (i in 1:ncol(file_paths))
{
  sample <- data.frame(read.table(file = file_paths[1,i], sep = '\t', header = TRUE))

  colnames(sample)[1] <- "sample_id"

  merge_table <- rbind(merge_table, sample)
}
# for loop that goes through the list of file paths for all samples and loads them
# then removes unnecessary and redundant columns then merges the sample's table by rows (samples)

merge_table <- drop_na(merge_table)

# rownames(merge_table) <- merge_table[1]
# causes error with incorrect number of rows (even though there is enough)

# merge_table <- merge_table[-1]

write.table(merge_table, file = snakemake@output[[1]], row.names=TRUE, sep="\t")
# save the merged data frame
