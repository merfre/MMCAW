### Script for creating species heatmap ###

library(ggplot2)
library("RColorBrewer")
library("gplots")

## Load files

sppno = snakemake.params.prevalence

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

blast_results <- blast_results[,order(colnames(blast_results))]
# order the columns by sample names for colors on plot

### If statement for determining most prevalent species

if(nrow(blast_results) <= sppno)
{
  ## Heatmap
  
  blast_abun <- data.frame(t(blast_results[order(rownames(blast_results)),]))
  
  lwid = c(1,3)
  lhei = c(1,3,1)
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  heatmap.2(t(data.matrix(blast_abun)), col = brewer.pal(n = 9, name = "Blues"),
            main = "All species identified by BLAST",
            margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
            key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
            trace="none", Rowv=FALSE, Colv=FALSE,
            cexCol = 1, cexRow = 1, srtCol = 45)
  
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
  
  ## Heatmap
  
  blast_top25_abun <- data.frame(t(blast_results[order(rownames(blast_results)),]))
  
  blast_top25_abun <- as.matrix(blast_top25_abun[colnames(blast_top25_abun) %in% blast_top25_names])
  
  lwid = c(1,3)
  lhei = c(1,3,1)
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  heatmap.2(t(data.matrix(blast_top25_abun)), col = brewer.pal(n = 9, name = "Blues"),
            main = paste("Top", sppno, "species identified by BLAST", sep=" "),
            margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
            key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
            trace="none", Rowv=FALSE, Colv=FALSE,
            cexCol = 1, cexRow = 1, srtCol = 45)
  
  dev.off()
  
}
