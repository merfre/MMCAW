### Script to count reads for each taxonomy of each sample ###

import pandas as pd

sample = pd.read_csv(snakemake.input[0], sep='\t')
# load sample results from LCA

path = snakemake.params.path
# load path
sample_no = path.split('/', -1) [1]
# separate sample name from path
sample = sample.drop(["query","classification","reason","lineage","lineage_scores"], axis=1).groupby(["species","genus","family","order","class","phylum","kingdom"]).size().to_frame(name = sample_no + "_cat").reset_index()
# drop unnecessary columns then group the dataframe by taxonomy and count the number of rows belonging
# to each taxonomy then turn it into a dataframe and make sample number the count column number
sample.to_csv(snakemake.output[0], sep='\t')
# save the final table
