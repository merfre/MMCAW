import pandas as pd
# version 2
# this uses merged.dmp to correct taxids from rankedlineage.dmp
# to use this, direct snakemake to the 'new_taxdump' directory rather than rankedlineage.dmp
# do this in the config file at TAXDUMP:

def reference_taxonomy():
    taxonomyDict = {}
    with open(snakemake.input.rankedlineage, 'r') as rankedlineage:
        for tax in rankedlineage:
            tax = tax.split('|')
            taxonid = tax[0]
            species = tax[1].strip().replace(
                ' ', '_') if tax[1].strip() else 'unidentified'
            genus = tax[3].strip() if tax[3].strip() else 'unidentified'
            family = tax[4].strip() if tax[4].strip() else 'unidentified'
            order = tax[5].strip() if tax[5].strip() else 'unidentified'
            classe = tax[6].strip() if tax[6].strip() else 'unidentified'
            phylum = tax[7].strip() if tax[7].strip() else 'unidentified'
            kingdom = tax[8].strip() if tax[8].strip() else 'unidentified'
            superkingdom = tax[9].strip() if tax[9].strip() else 'unidentified'
            taxonomyDict[str(tax[0].strip())] = {'species': species, 'genus': genus, 'family': family,
                                                 'order': order, 'class': classe, 'phylum': phylum, 'kingdom': kingdom, 'superkingdom': superkingdom}
    return taxonomyDict

def merged_taxonomy():
    mergedDict = {}
    with open(snakemake.input.merged, 'r') as merged:
        for taxid in merged:
            taxid = taxid.split('|')
            mergedDict[taxid[0].strip()] = taxid[1].strip()
    return mergedDict

taxonomyDict = reference_taxonomy()
mergedDict = merged_taxonomy()

with open(snakemake.input.kraken, 'r') as blasthits, open(snakemake.output.kraken_tax, 'w') as output:
    for hit in blasthits:
        taxid = hit.split('\t')[2].split(';')[0]
        if taxid == 'N/A' or taxid == '0':
            output.write(hit.strip() + '\tunidentified/unidentified/unidentified/unidentified/unidentified/unidentified/unidentified\n')
        else:
            try:
                superkingdom = taxonomyDict[taxid]['superkingdom']
            except KeyError:
                taxid = mergedDict[taxid]
                superkingdom = taxonomyDict[taxid]['superkingdom']

            if superkingdom != 'unidentified':
                output.write(hit.strip() + '\t' + taxonomyDict[taxid]['superkingdom']
                             + '/' + taxonomyDict[taxid]['phylum']
                             + '/' + taxonomyDict[taxid]['class']
                             + '/' + taxonomyDict[taxid]['order']
                             + '/' + taxonomyDict[taxid]['family']
                             + '/' + taxonomyDict[taxid]['genus']
                             + '/' + '_'.join(taxonomyDict[taxid]['species'].split('_')[:2]) + '\n')
            else:
                output.write(hit.strip() + '\t' + taxonomyDict[taxid]['kingdom']
                             + '/' + taxonomyDict[taxid]['phylum']
                             + '/' + taxonomyDict[taxid]['class']
                             + '/' + taxonomyDict[taxid]['order']
                             + '/' + taxonomyDict[taxid]['family']
                             + '/' + taxonomyDict[taxid]['genus']
                             + '/' + '_'.join(taxonomyDict[taxid]['species'].split('_')[:2]) + '\n')

# Need a step to separate the tax levels

sample = pd.read_csv(snakemake.output.kraken_tax, sep='\t')
# load sample results from LCA

sample.columns =['classified', 'query', 'taxid', 'length', 'score', 'taxa']
# add column names

# separate tax levels
taxa = sample['taxa'].str.split("/", n = -1, expand = True)
sample['kingdom'] = taxa[0]
sample['phylum'] = taxa[1]
sample['class'] = taxa[2]
sample['order'] = taxa[3]
sample['family'] = taxa[4]
sample['genus'] = taxa[5]
sample['species'] = taxa[6]
sample.drop(columns =["taxa"], inplace = True)

sample.to_csv(snakemake.output.kraken_tax, sep='\t')
# save the final table
