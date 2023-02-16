### R script to separate the taxonomy levels for the BLAST results merged in R ###

library(dplyr)
# Required for selecting certain columns

blast_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE)
# load the blast results for all samples across all taxonomies

### Defining function

separate_taxonomy <- function(out_path, tax_level) {

  blast_tax_level <- select(blast_results, ends_with("blast"), all_of(tax_level))
  # Separate samples and species levels

  tax_level_pos <- ncol(blast_tax_level)
  # Determine the column number of the taxonomy level

  write.table(blast_tax_level, file = out_path, row.names=TRUE, sep="\t")
}

### Species

separate_taxonomy(snakemake@output[[1]], "species")

### Family

#separate_taxonomy(snakemake@output[[2]], "family")

### Order

#separate_taxonomy(snakemake@output[[3]], "order")

### Class

#separate_taxonomy(snakemake@output[[4]], "class")

### Phyla

#separate_taxonomy(snakemake@output[[5]], "phyla")

### Kingdom

#separate_taxonomy(snakemake@output[[6]], "kingdom")
