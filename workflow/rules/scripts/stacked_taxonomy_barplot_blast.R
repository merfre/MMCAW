### Script to create side by side stacked barplots of taxonomy ###

library(ggplot2)
# for creating plots
library(tibble)
# for `rownames_to_column` and `column_to_rownames`
library(RColorBrewer)
# for plot colors

## Load files

sppno = snakemake@config[["prevalence"]]

# BLAST

blast_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

blast_results <- blast_results[!(blast_results$species=="unidentified" | blast_results$species=="unknown bacterium"),]
# remove results completely unclassified

row.names(blast_results) <- make.names(blast_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(blast_results)=="species"))
# determine which redundant columns to remove

blast_results <- blast_results[-columns_remove]
# Removes redundant species column

colnames(blast_results) <- gsub('_blast',' ',colnames(blast_results))
# Remove "_blast" suffix in column names

blast_results <- blast_results[,order(colnames(blast_results))]
# order the columns by sample names for colors on plot

### If statement for determining most prevalent species

if(nrow(blast_results) <= sppno)
{
  ## Barplot
  
  blast_abun <- data.frame(t(blast_results[order(rownames(blast_results)),]))
  
  sample_size <- nrow(blast_abun)
  species_color <- ncol(blast_abun)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  color=sample(col_vector, species_color)
  
  all_plot <- data.frame(
    Sample_ID=rep(c(rownames(blast_abun)), each = ncol(blast_abun)),
    Species=rep(c(colnames(blast_abun)), sample_size),
    Abundance=unlist(list(as.numeric(t(blast_abun))))
  )
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
    geom_bar(position="fill", stat="identity", colour="black")+
    scale_fill_manual(values = color)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    labs(title="BLAST taxonomy plot of all species identified"))
  #print(ggplot()) is required to prevent corrupted pdf in a loop/function
  dev.off()
  
} else
{
  ## Determine most prevalent species
  
  # BLAST
  
  logic_prop <- as.matrix(t(blast_results)[, 1:ncol(t(blast_results))] != 0)
  
  logic_prop <- data.frame(1*logic_prop)
  
  logic_prop <- data.frame(1*logic_prop)
  
  species_prop <- data.frame(colSums(logic_prop))
  
  species_propt <- data.frame(t(species_prop))
  
  sortspecies_propt <- sort(species_propt, decreasing=TRUE)
  
  blast_top25_species <- data.frame(sortspecies_propt[1:sppno])
  
  blast_top25_names <- colnames(blast_top25_species)
  
  blast_species_names <- colnames(sortspecies_propt)
  
  # Barplot
  
  blast_top25_abun <- data.frame(t(blast_results[order(rownames(blast_results)),]))
  
  blast_top25_abun <- as.matrix(blast_top25_abun[colnames(blast_top25_abun) %in% blast_top25_names])
  
  sample_size <- nrow(blast_top25_abun)
  species_color <- ncol(blast_top25_abun)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  color=sample(col_vector, species_color)
  
  all_plot <- data.frame(
    Sample_ID=rep(c(rownames(blast_top25_abun)), each = ncol(blast_top25_abun)),
    Species=rep(c(colnames(blast_top25_abun)), sample_size),
    Abundance=unlist(list(as.numeric(t(blast_top25_abun))))
  )
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
    geom_bar(position="fill", stat="identity", colour="black")+
    scale_fill_manual(values = color)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    labs(title=paste("BLAST taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
  #print(ggplot()) is required to prevent corrupted pdf in a loop/function
  dev.off()
  
}