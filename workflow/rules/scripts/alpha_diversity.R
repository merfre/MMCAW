### Script for alpha diversity analyses of CAT and Kraken2 results ###

library(vegan)
# package required for calculating alpha and beta diversities
library(dplyr)
# various commands to manipulate data frames
library(tibble)
# for `rownames_to_column` and `column_to_rownames`

## Load files

# Kraken2

kraken2_results <- read.table(file = snakemake@input[[2]], sep = '\t', header = TRUE, row.names = 1)

kraken2_results <- kraken2_results[!(kraken2_results$species=="unidentified"),]
# remove results completely unclassified

row.names(kraken2_results) <- make.names(kraken2_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(kraken2_results)=="species"))
# determine which redundant columns to remove

kraken2_results <- kraken2_results[-columns_remove]
# Removes redundant species column

if(ncol(kraken2_results) > 1)
{
  kraken2_results <- kraken2_results[,order(colnames(kraken2_results))]
  # order the columns by sample names for colors on plot
} else
{
  kraken2_results <- kraken2_results
}

# CAT

cat_results <- read.table(file = snakemake@input[[1]], sep = '\t', header = TRUE, row.names = 1)

cat_results <- cat_results[!(cat_results$species=="unidentified"),]
# remove results completely unclassified

row.names(cat_results) <- make.names(cat_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(cat_results)=="species"))
# determine which redundant columns to remove

cat_results <- cat_results[-columns_remove]
# Removes redundant species column

if(ncol(cat_results) > 1)
{
  cat_results <- cat_results[,order(colnames(cat_results))]
  # order the columns by sample names for colors on plot
} else
{
  cat_results <- cat_results
}

# BLAST

blast_results <- read.table(file = snakemake@input[[3]], sep = '\t', header = TRUE, row.names = 1)

blast_results <- blast_results[!(blast_results$species=="unidentified" | blast_results$species=="unknown bacterium"),]
# remove results completely unclassified

row.names(blast_results) <- make.names(blast_results$species, unique=TRUE)
# Changes row names to species names

columns_remove <- c(which(colnames(blast_results)=="species"))
# determine which redundant columns to remove

blast_results <- blast_results[-columns_remove]
# Removes redundant species column

if(ncol(blast_results) > 1)
{
  blast_results <- blast_results[,order(colnames(blast_results))]
  # order the columns by sample names for colors on plot
} else
{
  blast_results <- blast_results
}

## Alpha diversity table and plot

# CAT

invsimpson_cat <- diversity(t(cat_results),index = "invsimpson")
invsimpson_cat <- data.frame(t(invsimpson_cat))
row.names(invsimpson_cat) <- c("invsimpson_cat")
colnames(invsimpson_cat) <- colnames(cat_results)

evenness_cat <- eventstar(t(cat_results))
evenness_cat <- data.frame(t(evenness_cat))

shannon_cat <- diversity(t(cat_results), index="shannon")
shannon_cat <- data.frame(t(shannon_cat))
row.names(shannon_cat) <- c("shannon_cat")
colnames(shannon_cat) <- colnames(cat_results)

alpha_div_cat <- rbind(invsimpson_cat,evenness_cat[rownames(evenness_cat) == "Estar",],shannon_cat)
rownames(alpha_div_cat)[rownames(alpha_div_cat) == "2"] <- "evenness_cat"
colnames(alpha_div_cat) <- gsub('_cat','',colnames(alpha_div_cat))
alpha_div_cat <- rownames_to_column(alpha_div_cat, 'analysis')

# Kraken2

invsimpson_kraken <- diversity(t(kraken2_results),index = "invsimpson")
invsimpson_kraken <- data.frame(t(invsimpson_kraken))
row.names(invsimpson_kraken) <- c("invsimpson_kraken")
colnames(invsimpson_kraken) <- colnames(kraken2_results)

evenness_kraken <- eventstar(t(kraken2_results))
evenness_kraken <- data.frame(t(evenness_kraken))

shannon_kraken <- diversity(t(kraken2_results), index="shannon")
shannon_kraken <- data.frame(t(shannon_kraken))
row.names(shannon_kraken) <- c("shannon_kraken")
colnames(shannon_kraken) <- colnames(kraken2_results)

alpha_div_kraken <- rbind(invsimpson_kraken,evenness_kraken[rownames(evenness_kraken) == "Estar",],shannon_kraken)
rownames(alpha_div_kraken)[rownames(alpha_div_kraken) == "2"] <- "evenness_kraken"
colnames(alpha_div_kraken) <- gsub('_kraken','',colnames(alpha_div_kraken))
alpha_div_kraken <- rownames_to_column(alpha_div_kraken, 'analysis')

# BLAST

invsimpson_blast <- diversity(t(blast_results),index = "invsimpson")
invsimpson_blast <- data.frame(t(invsimpson_blast))
row.names(invsimpson_blast) <- c("invsimpson_blast")
colnames(invsimpson_blast) <- colnames(blast_results)

evenness_blast <- eventstar(t(blast_results))
evenness_blast <- data.frame(t(evenness_blast))

shannon_blast <- diversity(t(blast_results), index="shannon")
shannon_blast <- data.frame(t(shannon_blast))
row.names(shannon_blast) <- c("shannon_blast")
colnames(shannon_blast) <- colnames(blast_results)

alpha_div_blast <- rbind(invsimpson_blast,evenness_blast[rownames(evenness_blast) == "Estar",],shannon_blast)
rownames(alpha_div_blast)[rownames(alpha_div_blast) == "2"] <- "evenness_blast"
colnames(alpha_div_blast) <- gsub('_blast','',colnames(alpha_div_blast))
alpha_div_blast <- rownames_to_column(alpha_div_blast, 'analysis')

# Combine

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

fit <- lm(shannon_table[,2] ~ shannon_table[,1], data=data.frame(shannon_table))

pdf(snakemake@output[[2]])

plot(shannon_table[,c(1,2)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="CAT shannon score",
     main="Plot of shannon scores from BLAST and CAT results")
abline(fit, col = "black")

dev.off()

fit <- lm(shannon_table[,3] ~ shannon_table[,1], data=data.frame(shannon_table))

pdf(snakemake@output[[3]])

plot(shannon_table[,c(1,3)], type = "p", pch = 19, col = color, xlab="BLAST shannon score", ylab="Kraken2 shannon score",
     main="Plot of shannon scores from BLAST and Kraken2 results")
abline(fit, col = "black")

dev.off()

fit <- lm(shannon_table[,3] ~ shannon_table[,2], data=data.frame(shannon_table))

pdf(snakemake@output[[4]])

plot(shannon_table[,c(2,3)], type = "p", pch = 19, col = color, xlab="CAT shannon score", ylab="Kraken2 shannon score",
     main="Plot of shannon scores from CAT and Kraken2 results")
abline(fit, col = "black")

dev.off()

# Side by side shannon barplot

pdf(snakemake@output[[5]], width=25,height=20)

par(mar=c(11,4,4,4))
barplot(t(shannon_table), main="Barplot of BLAST, CAT, and Kraken2 shannon scores",
        xlab="Samples", ylab="Shannon Scores", col=c("blue","purple","red"),
        beside=TRUE, ylim=c(0,4), las=2)
legend(x=0,y=4, legend=c("BLAST","CAT","Kraken2"),
       fill=c("blue","purple","red"), bty="n", xpd=TRUE, cex=0.75, horiz=TRUE)

dev.off()
