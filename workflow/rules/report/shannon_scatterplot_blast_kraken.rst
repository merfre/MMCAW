Scatterplot of shannon scores from BLAST and Kraken2 results of all samples.



Assigner parameters:


BLAST was performed with a minimum percent identity of {{ snakemake.config["BLAST_min_perc_ident"] }},
minimum e-value of {{ snakemake.config["BLAST_min_evalue"] }}, and a max target sequence of {{ snakemake.config["BLAST_max_target_seqs"] }}.

Kraken2 was run with a confidence level of {{ snakemake.config["kraken_confidence"] }} and using the Kraken2 database "{{ snakemake.config["kraken_db"] }}"