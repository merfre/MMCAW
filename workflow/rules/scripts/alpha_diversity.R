### Script for alpha diversity analyses of CAT and Kraken2 results ###

library(vegan)
# package required for calculating alpha and beta diversities
library(dplyr)
# various commands to manipulate data frames
library(tibble)
# for `rownames_to_column` and `column_to_rownames`

## Transform data and create assigner alpha diversity table

alpha_div <- function(input, suffix) {
  
  ## Load files
  
  results <- read.table(file = input, sep = '\t', header = TRUE, row.names = 1)
  
  results <- results[!(results$species=="unidentified" | results$species=="unknown bacterium"),]
  # remove results completely unclassified
  
  row.names(results) <- make.names(results$species, unique=TRUE)
  # Changes row names to species names
  
  columns_remove <- c(which(colnames(results)=="species"))
  # determine which redundant columns to remove
  
  results <- results[-columns_remove]
  # Removes redundant species column
  
  if(ncol(results) > 1)
  {
    results <- results[,order(colnames(results))]
    # order the columns by sample names for colors on plot
  } else
  {
    results <- results
  }
  
  ## Alpha diversity table
  
  invsimpson <- diversity(t(results),index = "invsimpson")
  invsimpson <- data.frame(t(invsimpson))
  row.names(invsimpson) <- c(paste("invsimpson",suffix,sep=""))
  colnames(invsimpson) <- colnames(results)
  
  evenness <- eventstar(t(results))
  evenness <- data.frame(t(evenness))
  
  shannon <- diversity(t(results), index="shannon")
  shannon <- data.frame(t(shannon))
  row.names(shannon) <- c(paste("shannon",suffix,sep=""))
  colnames(shannon) <- colnames(results)
  
  output <- rbind(invsimpson,evenness[rownames(evenness) == "Estar",],shannon)
  rownames(output)[rownames(output) == "2"] <- paste("evenness",suffix,sep="")
  colnames(output) <- gsub(suffix,'',colnames(output))
  output <- rownames_to_column(output, 'analysis')
  
  return(output)
  
}

alpha_div_blast <- alpha_div(snakemake@input[[3]], '_blast')

alpha_div_cat <- alpha_div(snakemake@input[[1]], '_cat')

alpha_div_kraken <- alpha_div(snakemake@input[[2]], '_kraken')

## Alpha diversity plots

# Combine tables

alpha_div_table <- full_join(alpha_div_blast, alpha_div_cat)
alpha_div_table <- full_join(alpha_div_table, alpha_div_kraken)
alpha_div_table <- column_to_rownames(alpha_div_table, 'analysis')
# Put in alphabetical order (matters for scatterplots)

write.table(alpha_div_table, file = snakemake@output[[1]], row.names=TRUE, sep="\t")

# Shannon scatter plot

cat_rows <- grep("_cat", rownames(alpha_div_table))
kraken_rows <- grep("_kraken", rownames(alpha_div_table))
blast_rows <- grep("_blast", rownames(alpha_div_table))

shannon_rows <- grep("shannon", rownames(alpha_div_table))

shannon_table <- t(alpha_div_table[shannon_rows,])

sample_size <- ncol(alpha_div_table)

color <- c(replicate(sample_size, "red"))

# If statement depending on how many samples there are

if(sample_size > 1)
{
  
  # BLAST and CAT
  
  fit <- lm(shannon_table[,2] ~ shannon_table[,1], data=data.frame(shannon_table))
  
  pdf(snakemake@output[[2]])
  
  plot(shannon_table[,c(1,2)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="CAT shannon score",
       main="Plot of shannon scores from BLAST and CAT results")
  abline(fit, col = "black")
  
  dev.off()
  
  # Kraken2 and BLAST
  
  fit <- lm(shannon_table[,3] ~ shannon_table[,1], data=data.frame(shannon_table))
  
  pdf(snakemake@output[[3]])
  
  plot(shannon_table[,c(1,3)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="Kraken2 shannon score",
       main="Plot of shannon scores from BLAST and Kraken2 results")
  abline(fit, col = "black")
  
  dev.off()
  
  # CAT and Kraken2
  
  fit <- lm(shannon_table[,3] ~ shannon_table[,2], data=data.frame(shannon_table))
  
  pdf(snakemake@output[[4]])
  
  plot(shannon_table[,c(2,3)], type = "p", pch = 19, col = color, xlab="CAT shannon score", ylab="Kraken2 shannon score",
       main="Plot of shannon scores from CAT and Kraken2 results")
  abline(fit, col = "black")
  
  dev.off()

} else 
{
  
  # BLAST and CAT
  
  pdf(snakemake@output[[2]])
  
  plot(shannon_table[,c(1,2)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="CAT shannon score",
       main="Plot of shannon scores from BLAST and CAT results")
  
  dev.off()
  
  # Kraken2 and BLAST
  
  pdf(snakemake@output[[3]])
  
  plot(shannon_table[,c(1,3)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="Kraken2 shannon score",
       main="Plot of shannon scores from BLAST and Kraken2 results")
  
  dev.off()
  
  # CAT and Kraken2
  
  pdf(snakemake@output[[4]])
  
  plot(shannon_table[,c(2,3)], type = "p", pch = 19, col = color, xlab="CAT shannon score", ylab="Kraken2 shannon score",
       main="Plot of shannon scores from CAT and Kraken2 results")
  
  dev.off()
  
}

# Side by side shannon barplot

  pdf(snakemake@output[[5]], width=25,height=20)

  par(mar=c(11,4,4,4))
  barplot(t(shannon_table), main="Barplot of BLAST, CAT, and Kraken2 shannon scores",
          xlab="Samples", ylab="Shannon Scores", col=c("blue","purple","red"),
          beside=TRUE, ylim=c(0,4), las=2)
  legend(x=1,y=4, legend=c("BLAST","CAT","Kraken2"),
         fill=c("blue","purple","red"), bty="n", xpd=TRUE, cex=0.75, horiz=TRUE)

  dev.off()
