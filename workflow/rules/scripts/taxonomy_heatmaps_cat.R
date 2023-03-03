### Script for creating species heatmap ###

library(ggplot2)
library("RColorBrewer")
library("gplots")

## Load files

sppno = snakemake@config[["prevalence"]]

# cat

cat_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

cat_results <- cat_results[!(cat_results$species=="unidentified" | cat_results$species=="unknown bacterium"),]
# remove results completely unclassified

row.names(cat_results) <- make.names(cat_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(cat_results)=="species"))
# determine which redundant columns to remove

cat_results <- cat_results[-columns_remove]
# Removes redundant species column

colnames(cat_results) <- gsub('_cat',' ',colnames(cat_results))
# Remove "_cat" suffix in column names

if(ncol(cat_results) > 1)
{
  cat_results <- cat_results[,order(colnames(cat_results))]
  # order the columns by sample names for colors on plot
} else
{
  cat_results <- cat_results
}

### If statement for determining most prevalent species

if(nrow(cat_results) <= sppno)
{
  if(ncol(cat_results) > 1)
  {
    ## Heatmap
    
    cat_abun <- data.frame(t(cat_results[order(rownames(cat_results)),]))
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    heatmap.2(t(data.matrix(cat_abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = "Species identified by CAT",
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
    
  } else
  
  ## Heatmap
  
  cat_abun <- cat_results
  
  colnames(cat_abun) <- ""
  
  lwid = c(1,3)
  lhei = c(1,3,1)
  lmat = rbind(c(0,3),c(2,1),c(0,4))
  # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
  
  pdf(snakemake@output[[1]], width=25,height=20)
  
  heatmap.2(t(data.matrix(cbind(cat_abun[1], cat_abun[1]))), col = brewer.pal(n = 9, name = "Blues"),
            main = "Species identified by CAT",
            margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei, labRow = FALSE,
            key=T, key.title = "",key.ylab = "", key.xlab = "Read counts", ylab=colnames(cat_results),
            trace="none", Rowv=FALSE, Colv=FALSE,
            cexCol = 1, cexRow = 1, srtCol = 45)
  
  dev.off()
  
} else
{
  
  if(ncol(cat_results) > 1)
  {
    ## Determine most prevalent species
    
    # CAT
    
    logic_prop <- as.matrix(t(blast_results)[, 1:ncol(t(blast_results))] != 0)
    
    logic_prop <- data.frame(1*logic_prop)
    
    logic_prop <- data.frame(1*logic_prop)
    
    species_prop <- data.frame(colSums(logic_prop))
    
    species_propt <- data.frame(t(species_prop))
    
    sortspecies_propt <- sort(species_propt, decreasing=TRUE)
    
    cat_top25_species <- data.frame(sortspecies_propt[1:sppno])
    
    cat_top25_names <- colnames(cat_top25_species)
    
    cat_species_names <- colnames(sortspecies_propt)
    
    ## Heatmap
    
    cat_top25_abun <- data.frame(t(cat_results[order(rownames(cat_results)),]))
    
    cat_top25_abun <- as.matrix(cat_top25_abun[colnames(cat_top25_abun) %in% cat_top25_names])
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    heatmap.2(t(data.matrix(cat_top25_abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Top", sppno, "species identified by CAT", sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  } else
  {
    ## Determine most prevalent species
    
    # CAT
    
    logic_prop <- as.matrix(t(blast_results)[, 1:ncol(t(blast_results))] != 0)
    
    logic_prop <- data.frame(1*logic_prop)
    
    logic_prop <- data.frame(1*logic_prop)
    
    species_prop <- data.frame(colSums(logic_prop))
    
    species_propt <- data.frame(t(species_prop))
    
    sortspecies_propt <- sort(species_propt, decreasing=TRUE)
    
    cat_top25_species <- data.frame(sortspecies_propt[1:sppno])
    
    cat_top25_names <- colnames(cat_top25_species)
    
    cat_species_names <- colnames(sortspecies_propt)
    
    ## Heatmap
    
    cat_top25_abun <- data.frame(t(cat_results))
    
    cat_top25_abun <- as.matrix(cat_top25_abun[colnames(cat_top25_abun) %in% cat_top25_names])
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    heatmap.2(t(data.matrix(cat_top25_abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Top", sppno, "species identified by CAT", sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei, labCol = FALSE,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts", xlab=rownames(cat_top25_abun),
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  }
}