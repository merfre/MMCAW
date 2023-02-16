Table containing the taxonomy assignment comparison for {{ snakemake.wildcards["PATHS"] }} with each row being a sequence in that sample.

Contains all taxonomy levels and for each read the result of whether at least 2 taxonomy assignment tools were in agreement.
If no tool was able to assign a taxonomy the result is labelled 'unassigned' and if Kraken2, CAT, and BLAST were all in disagreement the result is labelled
'no_agreement'. Final 3 columns reports the percentage of taxonomies that each tool agreed on for that read.


Summary for all taxonomy comparison results can be found in the 'tax_comparison_summary' table.
