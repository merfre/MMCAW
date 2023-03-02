### Script to create side by side stacked barplots of taxonomy ###

library(ggplot2)
# for creating plots
library(tibble)
# for `rownames_to_column` and `column_to_rownames`
library(RColorBrewer)
# for plot colors

## Load files

sppno = snakemake@config[["prevalence"]]

# kraken

kraken_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

kraken_results <- kraken_results[!(kraken_results$species=="unidentified" | kraken_results$species=="unknown bacterium"),]
# remove results completely unclassified

row.names(kraken_results) <- make.names(kraken_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(kraken_results)=="species"))
# determine which redundant columns to remove

kraken_results <- kraken_results[-columns_remove]
# Removes redundant species column

colnames(kraken_results) <- gsub('_kraken',' ',colnames(kraken_results))
# Remove "_kraken" suffix in column names

kraken_results <- kraken_results[,order(colnames(kraken_results))]
# order the columns by sample names for colors on plot

### If statement for determining most prevalent species

if(nrow(kraken_results) <= sppno)
{
  ## Barplot
  
  kraken_abun <- data.frame(t(kraken_results[order(rownames(kraken_results)),]))
  
  sample_size <- nrow(kraken_abun)
  species_color <- ncol(kraken_abun)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  color=sample(col_vector, species_color)
  
  all_plot <- data.frame(
    Sample_ID=rep(c(rownames(kraken_abun)), each = ncol(kraken_abun)),
    Species=rep(c(colnames(kraken_abun)), sample_size),
    Abundance=unlist(list(as.numeric(t(kraken_abun))))
  )
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
    geom_bar(position="fill", stat="identity", colour="black")+
    scale_fill_manual(values = color)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    labs(title="Kraken2 taxonomy plot of all species identified"))
  #print(ggplot()) is required to prevent corrupted pdf in a loop/function
  dev.off()
  
} else
{
  ## Determine most prevalent species
  
  # kraken

  logic_prop <- as.matrix(t(kraken_results)[, 1:ncol(t(kraken_results))] != 0)

  logic_prop <- data.frame(1*logic_prop)

  logic_prop <- data.frame(1*logic_prop)

  species_prop <- data.frame(colSums(logic_prop))

  species_propt <- data.frame(t(species_prop))

  sortspecies_propt <- sort(species_propt, decreasing=TRUE)

  kraken_top25_species <- data.frame(sortspecies_propt[1:25])

  kraken_top25_names <- colnames(kraken_top25_species)

  kraken_species_names <- colnames(sortspecies_propt)

  # kraken

  kraken_top25_abun <- data.frame(t(kraken_results[order(rownames(kraken_results)),]))

  kraken_top25_abun <- as.matrix(kraken_top25_abun[colnames(kraken_top25_abun) %in% kraken_top25_names])

  sample_size <- nrow(kraken_top25_abun)
  species_color <- ncol(kraken_top25_abun)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  color=sample(col_vector, species_color)

  all_plot <- data.frame(
    Sample_ID=rep(c(rownames(kraken_top25_abun)), each = ncol(kraken_top25_abun)),
    Species=rep(c(colnames(kraken_top25_abun)), sample_size),
    Abundance=unlist(list(as.numeric(t(kraken_top25_abun))))
  )
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
    geom_bar(position="fill", stat="identity", colour="black")+
    scale_fill_manual(values = color)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    labs(title=paste("Kraken2 taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
  #print(ggplot()) is required to prevent corrupted pdf in a loop/function
  dev.off()
  
}