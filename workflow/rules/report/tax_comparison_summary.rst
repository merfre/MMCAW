Table containing summary information regarding direct comparison of taxonomy assignment results from Kraken2, CAT, and BLAST.

Each row is a sample and contains the percentage of taxonomy assignments that 2 or more tools agreed on.
Also contains the overall percentage that each tool agreed with eachother on taxonomy assignments for that sample. 



Assigner parameters:


CAT was run using default parameters, the CAT database "{{ snakemake.config["cat_db"] }}",
and the CAT taxonomy reference "{{ snakemake.config["cat_taxonomy"] }}".

BLAST was performed with a minimum percent identity of {{ snakemake.config["BLAST_min_perc_ident"] }},
minimum e-value of {{ snakemake.config["BLAST_min_evalue"] }}, and a max target sequence of {{ snakemake.config["BLAST_max_target_seqs"] }}.

Kraken2 was run with a confidence level of {{ snakemake.config["kraken_confidence"] }} and using the Kraken2 database "{{ snakemake.config["kraken_db"] }}"