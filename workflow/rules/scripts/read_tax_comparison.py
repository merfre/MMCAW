### Script to compare the taxonomies from all assigners for each read ###

import pandas as pd
import numpy as np
# load packages

kraken = pd.read_csv(snakemake.input.kraken_results, sep='\t').replace("unknown", "unidentified")
# load kraken
cat = pd.read_csv(snakemake.input.cat_results, sep='\t').replace("unknown", "unidentified")
# load cat
blast = pd.read_csv(snakemake.input.blast_results, sep='\t').replace("unknown", "unidentified")
# load blast
# Looking at blast output I'm seeing unknown in the outputs, so all "unknown" is replaced with "unidentified" to maintain uniformity

# separate query and tax levels for all 3 inputs
kraken = kraken[['query','kingdom','phylum','class','order','family','genus','species',]]
cat = cat[['query','kingdom','phylum','class','order','family','genus','species',]]
blast = blast[['query','kingdom','phylum','class','order','family','genus','species',]]

# so that no reads are missed, a list of all query names is named so that assigner results that did not recieve hits for a read
# can include the missing read in their results
kraken_query = list(kraken['query'])
cat_query = list(cat['query'])
blast_query = list(blast['query'])
query_list = pd.DataFrame(list(set(kraken_query + cat_query + blast_query)))
query_list.columns = ['query']

kraken = query_list.merge(kraken, how='left',on='query').replace(np.nan, 'unassigned', regex=True)
cat = query_list.merge(cat, how='left',on='query').replace(np.nan, 'unassigned', regex=True)
blast = query_list.merge(blast, how='left',on='query').replace(np.nan, 'unassigned', regex=True)

# query is made the index (row names) so that the dataframes can be compared by each row
# then the redundunt query column is dropped

# kraken
kraken.index = kraken['query']
kraken = kraken.drop(['query',], axis=1).sort_index(ascending=True)
# cat
cat.index = cat['query']
cat = cat.drop(['query',], axis=1).sort_index(ascending=True)
# blast
blast.index = blast['query']
blast = blast.drop(['query',], axis=1).sort_index(ascending=True)

# After examining the outputs and NCBI taxonomies it's clear that when Kraken and CAT cannot find an assignment they label
# the read unidentified at every level (kingdom is the real marker for a read with no assignment)

# this function looks for whether all the columns in each row are equal, which only occurs when the read is unidentified at
# every level and therefore unassigned
kraken_equal = kraken.eq(kraken.iloc[:, 0], axis=0).all(1)
cat_equal = cat.eq(cat.iloc[:, 0], axis=0).all(1)

def unassigned(arg1, equality):
    for read in range(len(equality)) :
        if equality.iloc[read]:
            arg1.iloc[read, :] = 'unassigned'
        else:
            arg1.iloc[read, :] = arg1.iloc[read, :]

# The function goes through the list of True/False for the assigner (whether the columns of the row contain the same value)
# and if the read has True (all columns the same) then all the columns of that row are replaced with unassigned
# if the read has false in the equality list then all the columns are left the same

unassigned(kraken, kraken_equal)
unassigned(cat, cat_equal)

# create empty results tables with appropriate index and column names

column_names = ['kingdom','phylum','class','order','family','genus','species',]

results_table = pd.DataFrame(columns = column_names, index = query_list['query']).sort_index(ascending=True)
results_table2 = pd.DataFrame(columns = column_names, index = query_list['query']).sort_index(ascending=True)

# Function to determine how much all 3 outputs agree

def agree3(arg1, arg2, arg3, level, output):
    colno = arg1.columns.get_loc(level)
    for read in range(len(arg1)) :
        if arg1.iloc[read, colno] == 'unassigned':
            if arg2.iloc[read, colno] == arg3.iloc[read, colno]:
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg3.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg2.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg3.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg2.iloc[read, colno] == 'unassigned':
            if arg1.iloc[read, colno] == arg3.iloc[read, colno]:
                output.iloc[read, colno] = arg1.iloc[read, colno]
            elif arg3.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg1.iloc[read, colno]
            elif arg1.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg3.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg3.iloc[read, colno] == 'unassigned':
            if arg2.iloc[read, colno] == arg1.iloc[read, colno]:
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg1.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg2.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg1.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg1.iloc[read, colno] != arg2.iloc[read, colno]:
            output.iloc[read, colno] = 'no_agreement'
        elif arg2.iloc[read, colno] != arg3.iloc[read, colno]:
            output.iloc[read, colno] = 'no_agreement'
        elif arg1.iloc[read, colno] != arg3.iloc[read, colno]:
            output.iloc[read, colno] = 'no_agreement'
        else:
            output.iloc[read, colno] = arg1.iloc[read, colno] # therefore kraken and all outputs are equal
# first the column number of the level in question is saved as colno then a for loop to go through each read (each row of the dataframes)
# then each pair of assigners are compared and if they are not the same then the cell is filled with the kraken output since all results are the same
# the .iloc function and the use of numbers to locate the columns is necessary because one of the columns is called class
# if numbers are not used to locate each cell of the table, an error is thrown when accessing the class column

agree3(kraken, cat, blast, 'kingdom', results_table2)
agree3(kraken, cat, blast, 'phylum', results_table2)
agree3(kraken, cat, blast, 'class', results_table2)
agree3(kraken, cat, blast, 'order', results_table2)
agree3(kraken, cat, blast, 'family', results_table2)
agree3(kraken, cat, blast, 'genus', results_table2)
agree3(kraken, cat, blast, 'species', results_table2)
# results_table2 is mostly filled with 'no_agreement' and very few results are agreed upon with all 3 assigners

# Function to vote via majority rule for the assignment of the read without counting unassigned as a vote

def agree2(arg1, arg2, arg3, level, output):
    colno = arg1.columns.get_loc(level)
    for read in range(len(arg1)) :
        if arg1.iloc[read, colno] == arg2.iloc[read, colno] and arg2.iloc[read, colno] == arg3.iloc[read, colno] and arg1.iloc[read, colno] == arg3.iloc[read, colno] :
            output.iloc[read, colno] = arg1.iloc[read, colno]
        elif arg1.iloc[read, colno] == 'unassigned':
            if arg2.iloc[read, colno] == arg3.iloc[read, colno]:
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg3.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg2.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg3.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg2.iloc[read, colno] == 'unassigned':
            if arg1.iloc[read, colno] == arg3.iloc[read, colno]:
                output.iloc[read, colno] = arg1.iloc[read, colno]
            elif arg3.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg1.iloc[read, colno]
            elif arg1.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg3.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg3.iloc[read, colno] == 'unassigned':
            if arg2.iloc[read, colno] == arg1.iloc[read, colno]:
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg1.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg2.iloc[read, colno]
            elif arg2.iloc[read, colno] == 'unassigned':
                output.iloc[read, colno] = arg1.iloc[read, colno]
            else:
                output.iloc[read, colno] = 'no_agreement'
        elif arg1.iloc[read, colno] == arg2.iloc[read, colno]:
            output.iloc[read, colno] = arg1.iloc[read, colno]
        elif arg2.iloc[read, colno] == arg3.iloc[read, colno]:
            output.iloc[read, colno] = arg2.iloc[read, colno]
        elif arg1.iloc[read, colno] == arg3.iloc[read, colno]:
            output.iloc[read, colno] = arg1.iloc[read, colno]
        else:
            output.iloc[read, colno] = 'no_agreement'
# first the column number of the level in question is saved as colno then a for loop to go through each read (each row of the dataframes)
# Then if all 3 assigners agree then the results table cell is marked as argument's one result (since they're all the same)
# next a pairwise comparison is performed for all the assigners to ignore the unassigned result unless all 3 report unassigned
# next each assigner output is compared and if any agree then the result is saved
# the .iloc function and the use of numbers to locate the columns is necessary because one of the columns is called class
# if numbers are not used to locate each cell of the table, an error is thrown when accessing the class column

agree2(kraken, cat, blast, 'kingdom', results_table)
agree2(kraken, cat, blast, 'phylum', results_table)
agree2(kraken, cat, blast, 'class', results_table)
agree2(kraken, cat, blast, 'order', results_table)
agree2(kraken, cat, blast, 'family', results_table)
agree2(kraken, cat, blast, 'genus', results_table)
agree2(kraken, cat, blast, 'species', results_table)

# create an empty table for holding the individual agreement for all the assigners with appropriate index and column names

column_names = ['kingdom','phylum','class','order','family','genus','species',]

blast_cat_agreement = pd.DataFrame(columns = column_names, index = query_list['query']).sort_index(ascending=True)
kraken_blast_agreement = pd.DataFrame(columns = column_names, index = query_list['query']).sort_index(ascending=True)
cat_kraken_agreement = pd.DataFrame(columns = column_names, index = query_list['query']).sort_index(ascending=True)

# blast:cat, kraken:blast, cat:kraken

# Function to count how much each assigner agrees for each read

def agree_each(arg1, arg2, output):
    for read in range(len(arg1)) :
        if arg1.iloc[read,0] == arg2.iloc[read,0]:
            output.iloc[read,0] = 1
        else:
            output.iloc[read,0] = 0
        if arg1.iloc[read,2] == arg2.iloc[read,2]:
            output.iloc[read,2] = 1
        else:
            output.iloc[read,2] = 0
        if arg1.iloc[read,3] == arg2.iloc[read,3]:
            output.iloc[read,3] = 1
        else:
            output.iloc[read,3] = 0
        if arg1.iloc[read,4] == arg2.iloc[read,4]:
            output.iloc[read,4] = 1
        else:
            output.iloc[read,4] = 0
        if arg1.iloc[read,5] == arg2.iloc[read,5]:
            output.iloc[read,5] = 1
        else:
            output.iloc[read,5] = 0
        if arg1.iloc[read,6] == arg2.iloc[read,6]:
            output.iloc[read,6] = 1
        else:
            output.iloc[read,6] = 0
        if arg1.iloc[read,1] == arg2.iloc[read,1]:
            output.iloc[read,1] = 1
        else:
            output.iloc[read,1] = 0
# first there's a for loop to go through all the reads then for each read, at each taxonomy level the 2 assigner results are compared
# if they agree a one is saved, if they don't a 0 is saved in the blank table previously created
# the .iloc function and the use of numbers to locate the columns is necessary because one of the columns is called class
# if numbers are not used to locate each cell of the table, an error is thrown when accessing the class column

agree_each(blast,cat,blast_cat_agreement)
agree_each(kraken,blast,kraken_blast_agreement)
agree_each(cat,kraken,cat_kraken_agreement)

# Add blank column for each assigner comparison

results_table["blast_cat_agreement"] = ""
results_table["kraken_blast_agreement"] = ""
results_table["cat_kraken_agreement"] = ""

# for each query (read) get the average of tax levels that the 2 assigners agreed on
# then get the average number of levels the 2 assigners agreed on and multiple by 100 to get the percentage that the 2 assigners for that read

for read in range(len(blast_cat_agreement)) :
    results_table.iloc[read,7] = blast_cat_agreement.iloc[read,0:].mean() * 100
    results_table.iloc[read,8] = kraken_blast_agreement.iloc[read,0:].mean() * 100
    results_table.iloc[read,9] = cat_kraken_agreement.iloc[read,0:].mean() * 100

# create an empty table for summary stats with appropriate index (sample ID) and column names

sample_ID = snakemake.params.path.rsplit('/', -1)[1]

column_names = ['kingdom_no_agreements','phylum_no_agreements','class_no_agreements','order_no_agreements','family_no_agreements','genus_no_agreements','species_no_agreements','blast_cat_agreement','kraken_blast_agreement','cat_kraken_agreement',]

tax_comparison_summary = pd.DataFrame(columns = column_names, index=[sample_ID,])

# for loop to go through the first 6 columns of results and summary table and count how mant results are no_agreement

for col in range(7) :
        if results_table.iloc[:,col].str.find('no_agreement').sum() == -len(results_table):
            tax_comparison_summary.iloc[:,col] = 0
        else:
            tax_comparison_summary.iloc[:,col] = results_table.iloc[:,col].value_counts().no_agreement / len(results_table) * 100

# if statement is required because if the column contains no no_agreement results an error occurs
# to prevent the error if statement checks whether no_agreement is in the column, if no_agreement is missing it is
# reported as -1 ; therefore the total result would be negative number of rows and that result would be reported as 0
# If any column contains a no_agreement result the number of no_agreements is counted, divided by the number of reads,
# multiplied by 100 (to get percentage of no_agreements in the column), and saved in the blank table

tax_comparison_summary['blast_cat_agreement'] = results_table['blast_cat_agreement'].mean()
tax_comparison_summary['kraken_blast_agreement'] = results_table['kraken_blast_agreement'].mean()
tax_comparison_summary['cat_kraken_agreement'] = results_table['cat_kraken_agreement'].mean()

tax_comparison_summary.to_csv(snakemake.output.sample_summary,sep="\t")
results_table.to_csv(snakemake.output.sample_comparison,sep="\t")
