### Script for beta diversity analyses of CAT and Kraken2 results ###

library(vegan)
# package required for calculating alpha and beta diversities
library(stringi)
# required for ggfortify and ggplot2
library(RColorBrewer)
# package for colors for plots
library(ggplot2)
# package for making plots
library(ggfortify)
# package that ggplot2 requires for making PCA plot
library(gdata)
# To combine the taxonomy assigner results for distance matrices
library(ggrepel)
library(dplyr)

## Load files

# Kraken2

kraken2_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

kraken2_results <- kraken2_results[,order(colnames(kraken2_results))]
# order the columns by sample names for colors on plot

kraken2_results <- aggregate(.~species, kraken2_results,FUN=sum)
# sums rows of the same species name

# CAT

cat_results <- read.table(file = snakemake@input[[2]], sep = '\t', header = TRUE, row.names = 1)

cat_results <- cat_results[,order(colnames(cat_results))]
# order the columns by sample names for colors on plot

cat_results <- aggregate(.~species, cat_results,FUN=sum)
# sums rows of the same species name

# BLAST

blast_results <- read.table(file = snakemake@input[[3]], sep = '\t', header = TRUE, row.names = 1)

blast_results <- blast_results[,order(colnames(blast_results))]
# order the columns by sample names for colors on plot

blast_results <- aggregate(.~species, blast_results,FUN=sum)
# sums rows of the same species name

## Beta diversity PCA plot

combined <- full_join(kraken2_results, cat_results)
# full_join combines the results despite different number of rows

combined <- full_join(blast_results, combined)

row.names(combined) <- make.names(combined$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(combined)=="species"))
# determine which redundant columns to remove

combined <- combined[-columns_remove]
# Removes redundant species column

combined <- data.frame(t(combined))

combined[is.na(combined)] <- 0
# this step is required because not all methods found the same species and this keeps these species in the analysis

combined_vegd <- vegdist(na.omit(combined), method="bray", na.rm = TRUE)

combined_vegpca <- prcomp(combined_vegd, center=TRUE, scale=TRUE)

# Colors for samples

sample_size <- ncol(kraken2_results) - 1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

color=sample(col_vector, sample_size, replace = TRUE)

color_all <- data.frame(seq(c(2*sample_size)))

for(i in 1:sample_size)
{
  color_all[i,1] <- color[i]
}

c = sample_size

for(i in 1:sample_size)
{
  c = c + 1

  color_all[c,1] <- color[i]
}

c = sample_size*2

for(i in 1:sample_size)
{
  c = c + 1

  color_all[c,1] <- color[i]
}

# Making plots

pdf(snakemake@output[[1]], width=25,height=20)

autoplot(combined_vegpca, colour = as.matrix(color_all), label=TRUE, label.repel=TRUE,
         main ="PCA of species Beta diversity from results of all taxonomic assigners")

dev.off()
