### Script to combine and merge blast results into 1 table ###

library(dplyr)
# Package required to get distinct taxonomies and levels

file_paths <- t(read.csv(snakemake@input[[1]], header=FALSE))
# Load a list of file paths of the individual CAT results

taxonomy_list <- data.frame(read.table(file = snakemake@input[[2]], sep = '\t', header = FALSE))
# Load the list of all unique clades identified by cat in this batch of samples

colnames(taxonomy_list) <- c("species","genus","family","order","class","phylum","kingdom")
# Add column names to taxonomy list so that data frames can be merged based on column names

taxonomy_list <- taxonomy_list[!(taxonomy_list$species=="species"),]
# Remove the column names that were included in the list of taxonomies

new_col = ncol(taxonomy_list)

for (i in 1:ncol(file_paths))
{
  sample <- data.frame(read.table(file = file_paths[1,i], sep = '\t', header = TRUE))

  sample <- sample[,-1]

  taxonomy_list <- merge(taxonomy_list, sample, by=c("species","genus","family","order","class","phylum","kingdom"), all.x=T)

  taxonomy_list[[new_col + i]][is.na(taxonomy_list[[new_col + i]])] = 0
}
# for loop that goes through the list of file paths for all samples and loads them
# then removes unnecessary and redundant columns then merges the sample's table by all taxonomy levels
# finally the NAs in the new column are replaced with 0

write.table(taxonomy_list, file = snakemake@output[[1]], row.names=TRUE, sep="\t")
# save the merged data frame
