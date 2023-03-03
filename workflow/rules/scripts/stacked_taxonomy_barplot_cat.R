### Script to create side by side stacked barplots of taxonomy ###

library(ggplot2)
# for creating plots
library(tibble)
# for `rownames_to_column` and `column_to_rownames`
library(RColorBrewer)
# for plot colors

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
    ## Barplot
    
    cat_abun <- data.frame(t(cat_results[order(rownames(cat_results)),]))
    
    sample_size <- nrow(cat_abun)
    species_color <- ncol(cat_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(cat_abun)), each = ncol(cat_abun)),
      Species=rep(c(colnames(cat_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(cat_abun))))
    )
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title="CAT taxonomy plot of all species identified"))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  } else
  {
    ## Barplot
    
    cat_abun <- data.frame(t(cat_results))
    
    sample_size <- nrow(cat_abun)
    species_color <- ncol(cat_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(cat_abun)), each = ncol(cat_abun)),
      Species=rep(c(colnames(cat_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(cat_abun))))
    )
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title="CAT taxonomy plot of all species identified"))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  }
} else
{
  if(ncol(cat_results) > 1)
  {
    ## Determine most prevalent species
    
    # cat
    
    logic_prop <- as.matrix(t(cat_results)[, 1:ncol(t(cat_results))] != 0)
    
    logic_prop <- data.frame(1*logic_prop)
    
    logic_prop <- data.frame(1*logic_prop)
    
    species_prop <- data.frame(colSums(logic_prop))
    
    species_propt <- data.frame(t(species_prop))
    
    sortspecies_propt <- sort(species_propt, decreasing=TRUE)
    
    cat_top25_species <- data.frame(sortspecies_propt[1:sppno])
    
    cat_top25_names <- colnames(cat_top25_species)
    
    cat_species_names <- colnames(sortspecies_propt)
    
    # cat
    
    cat_top25_abun <- data.frame(t(cat_results[order(rownames(cat_results)),]))
    
    cat_top25_abun <- as.matrix(cat_top25_abun[colnames(cat_top25_abun) %in% cat_top25_names])
    
    sample_size <- nrow(cat_top25_abun)
    species_color <- ncol(cat_top25_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(cat_top25_abun)), each = ncol(cat_top25_abun)),
      Species=rep(c(colnames(cat_top25_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(cat_top25_abun))))
    )
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste("CAT taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  } else
  {
    ## Determine most prevalent species
    
    # cat
    
    logic_prop <- as.matrix(t(cat_results)[, 1:ncol(t(cat_results))] != 0)
    
    logic_prop <- data.frame(1*logic_prop)
    
    logic_prop <- data.frame(1*logic_prop)
    
    species_prop <- data.frame(colSums(logic_prop))
    
    species_propt <- data.frame(t(species_prop))
    
    sortspecies_propt <- sort(species_propt, decreasing=TRUE)
    
    cat_top25_species <- data.frame(sortspecies_propt[1:sppno])
    
    cat_top25_names <- colnames(cat_top25_species)
    
    cat_species_names <- colnames(sortspecies_propt)
    
    # cat
    
    cat_top25_abun <- data.frame(t(cat_results))
    
    cat_top25_abun <- as.matrix(cat_top25_abun[colnames(cat_top25_abun) %in% cat_top25_names])
    
    sample_size <- nrow(cat_top25_abun)
    species_color <- ncol(cat_top25_abun)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    color=sample(col_vector, species_color)
    
    all_plot <- data.frame(
      Sample_ID=rep(c(rownames(cat_top25_abun)), each = ncol(cat_top25_abun)),
      Species=rep(c(colnames(cat_top25_abun)), sample_size),
      Abundance=unlist(list(as.numeric(t(cat_top25_abun))))
    )
    
    pdf(snakemake@output[[1]], width=25,height=20)
    
    print(ggplot(all_plot, aes(fill=Species, y=Abundance, x=Sample_ID)) +
            geom_bar(position="fill", stat="identity", colour="black")+
            scale_fill_manual(values = color)+
            theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                             size = 12, hjust = 1))+
            labs(title=paste("CAT taxonomy plot of the top", sppno, "most prevalent species", sep=" ")))
    #print(ggplot()) is required to prevent corrupted pdf in a loop/function
    dev.off()
  }
}