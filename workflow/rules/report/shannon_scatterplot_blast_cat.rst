Scatterplot of shannon scores from BLAST and CAT results of all samples.



Assigner parameters:


CAT was run using default parameters, the CAT database "{{ snakemake.config["cat_db"] }}",
and the CAT taxonomy reference "{{ snakemake.config["cat_taxonomy"] }}".

BLAST was performed with a minimum percent identity of {{ snakemake.config["BLAST_min_perc_ident"] }},
minimum e-value of {{ snakemake.config["BLAST_min_evalue"] }}, and a max target sequence of {{ snakemake.config["BLAST_max_target_seqs"] }}.
