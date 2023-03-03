### Script for creating species heatmap ###

library(ggplot2)
library("RColorBrewer")
library("gplots")

## Load files

sppno = snakemake@config[["prevalence"]]

# cat

cat_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)


### Defining function

taxonomy_heatmap <- function(input, assigner, suffix, out_path) {
  
  input <- input[!(input$species=="unidentified" | input$species=="unknown bacterium"),]
  # remove results completely unclassified
  
  row.names(input) <- make.names(input$species, unique=TRUE)
  # Changes row names to species names
  
  columns_remove <- c(which(colnames(input)=="species"))
  # determine which redundant columns to remove
  
  input <- input[-columns_remove]
  # Removes redundant species column
  
  colnames(input) <- gsub(suffix,' ',colnames(input))
  # Remove assigner suffix in column names
  
  if(ncol(input) > 1)
  {
    input <- input[,order(colnames(input))]
    # order the columns by sample names for colors on plot
  } else
  {
    input <- input
  }
  
  ### If statements for determining number of samples and most prevalent species
  
  if(nrow(input) >= sppno & ncol(input) > 1)
  {
    ## Determine most prevalent species
    
    logic_prop <- as.matrix(t(input)[, 1:ncol(t(input))] != 0)
    
    logic_prop <- data.frame(1*logic_prop)
    
    logic_prop <- data.frame(1*logic_prop)
    
    species_prop <- data.frame(colSums(logic_prop))
    
    species_propt <- data.frame(t(species_prop))
    
    sortspecies_propt <- sort(species_propt, decreasing=TRUE)
    
    top_species <- data.frame(sortspecies_propt[1:sppno])
    
    top_names <- colnames(top_species)
    
    species_names <- colnames(sortspecies_propt)
    
    ## Heatmap
    
    top_abun <- data.frame(t(input[order(rownames(input)),]))
    
    top_abun <- as.matrix(top_abun[colnames(top_abun) %in% top_names])
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(out_path, width=25,height=20)
    
    heatmap.2(t(data.matrix(top_abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Top", sppno, "species identified by", assigner, sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  } else if (nrow(input) < sppno & ncol(input) > 1)
  {
    ## Heatmap
    
    abun <- data.frame(t(input[order(rownames(input)),]))
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(out_path, width=25,height=20)
    
    heatmap.2(t(data.matrix(abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Species identified by", assigner, sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts",
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  } else if (nrow(input) >= sppno & ncol(input) < 1)
  {
    ## Determine most abundant species
    
    sortspecies <- sort(data.frame(t(input)), decreasing=TRUE)
    
    top_species <- data.frame(sortspecies[1:sppno])
    
    top_names <- colnames(top_species)
    
    species_names <- colnames(sortspecies)
    
    ## Heatmap
    
    top_abun <- data.frame(t(input))
    
    top_abun <- as.matrix(top_abun[colnames(top_abun) %in% top_names])
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(out_path, width=25,height=20)
    
    heatmap.2(t(data.matrix(top_abun)), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Top", sppno, "species identified by", assigner, sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei, labCol = FALSE,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts", xlab=rownames(top_abun),
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  } else (nrow(input) < sppno & ncol(input) < 1)
  {
    ## Heatmap
    
    abun <- input
    
    colnames(abun) <- ""
    
    lwid = c(1,3)
    lhei = c(1,3,1)
    lmat = rbind(c(0,3),c(2,1),c(0,4))
    # Moves the key to the bottom of the figure and ensures the whole title and heatmap are included
    
    pdf(out_path, width=25,height=20)
    
    heatmap.2(t(data.matrix(cbind(abun[1], abun[1]))), col = brewer.pal(n = 9, name = "Blues"),
              main = paste("Species identified by", assigner, sep=" "),
              margins = c(10,13),lmat = lmat, lwid = lwid, lhei = lhei, labRow = FALSE,
              key=T, key.title = "",key.ylab = "", key.xlab = "Read counts", ylab=colnames(input),
              trace="none", Rowv=FALSE, Colv=FALSE,
              cexCol = 1, cexRow = 1, srtCol = 45)
    
    dev.off()
  }
}

### Plot creation

taxonomy_heatmap(cat_results, "CAT", '_cat', snakemake@output[[1]])
