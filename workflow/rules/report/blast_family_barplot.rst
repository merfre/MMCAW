Barplot of family level assignment from BLAST results.


BLAST was performed with a minimum percent identity of {{ snakemake.config["BLAST_min_perc_ident"] }},
minimum e-value of {{ snakemake.config["BLAST_min_evalue"] }}, and a max target sequence of {{ snakemake.config["BLAST_max_target_seqs"] }}.

The LCA assignment of BLAST results used the BLAST database "{{ snakemake.config["blast_db"] }}", a bitscore threshold of {{ snakemake.config["MLCA_bitscore"] }},
identity of {{ snakemake.config["MLCA_identity"] }}, coverage of {{ snakemake.config["MCLA_coverage"] }},
majority score of {{ snakemake.config["MLCA_majority"] }}, and a hit threshold of {{ snakemake.config["MLCA_hits"] }}.
