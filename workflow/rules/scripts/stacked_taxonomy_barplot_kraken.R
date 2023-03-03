### Script to create side by side stacked barplots of taxonomy ###

library(ggplot2)
# for creating plots
library(tibble)
# for `rownames_to_column` and `column_to_rownames`
library(RColorBrewer)
# for plot colors

## Load files

sppno = snakemake@config[["prevalence"]]

# Kraken2

kraken_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

### Defining function

taxonomy_stacked_barplot <- function(input, assigner, suffix, out_path) {
  
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
    
    ## Barplot
    
    top_abun <- data.frame(t(input[order(rownames(input)),]))
    
    top_abun <- as.matrix(top_abun[colnames(top_abun) %in% top_names])
    
    sample_size <- nrow(top_abun)
    species_color <- ncol(top_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(top_abun)), each = ncol(top_abun)),
      Species=rep(c(colnames(top_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(top_abun))))
    )
    
    pdf(out_path, width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste(assigner, "taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  } else if (nrow(input) < sppno & ncol(input) > 1)
  {
    ## Barplot
    
    abun <- data.frame(t(input[order(rownames(input)),]))
    
    sample_size <- nrow(abun)
    species_color <- ncol(abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(abun)), each = ncol(abun)),
      Species=rep(c(colnames(abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(abun))))
    )
    
    pdf(out_path, width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste(assigner, "taxonomy plot of all species identified", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  } else if (nrow(input) >= sppno & ncol(input) < 1)
  {
    ## Determine most abundant species
    
    sortspecies <- sort(data.frame(t(input)), decreasing=TRUE)
    
    top_species <- data.frame(sortspecies[1:sppno])
    
    top_names <- colnames(top_species)
    
    species_names <- colnames(sortspecies)
    
    ## Barplot
    
    top_abun <- data.frame(t(input))
    
    top_abun <- as.matrix(top_abun[colnames(top_abun) %in% top_names])
    
    sample_size <- nrow(top_abun)
    species_color <- ncol(top_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(top_abun)), each = ncol(top_abun)),
      Species=rep(c(colnames(top_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(top_abun))))
    )
    
    pdf(out_path, width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste(assigner, "taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  } else (nrow(input) < sppno & ncol(input) < 1)
  {
    ## Barplot
    
    abun <- data.frame(t(input))
    
    sample_size <- nrow(abun)
    species_color <- ncol(abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(abun)), each = ncol(abun)),
      Species=rep(c(colnames(abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(abun))))
    )
    
    pdf(out_path, width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste(assigner, "taxonomy plot of all species identified", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  }
}

### Plot creation

taxonomy_stacked_barplot(kraken_results, "Kraken2", '_kraken', snakemake@output[[1]])