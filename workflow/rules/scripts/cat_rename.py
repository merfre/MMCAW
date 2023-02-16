### Script to rename the columns in CAT summary output

import pandas as pd

sample = pd.read_csv(snakemake.input[0], sep='\t', header=None)

path = snakemake.params[0]
# load path
sample_no = path.split('/', -1) [1]
# separate sample name from path
sample = sample.replace(regex='number of contigs', value=sample_no)
# replace the column name "number of contigs" with the sample name so that the samples can be combined in R
sample.columns = sample.iloc[0]
# make the first column the header
sample = sample[1:]
# take the data less the header row
sample.to_csv(snakemake.output[0], sep='\t')
# save the final table
