Scatterplot of shannon scores from Kraken2 and CAT results of all samples.



Assigner parameters:


CAT was run using default parameters, the CAT database "{{ snakemake.config["cat_db"] }}",
and the CAT taxonomy reference "{{ snakemake.config["cat_taxonomy"] }}".

Kraken2 was run with a confidence level of {{ snakemake.config["kraken_confidence"] }} and using the Kraken2 database "{{ snakemake.config["kraken_db"] }}"